import pandas as pd
from typing import Dict
import vcfpy
import numpy as np
import os
import bioframe as bf
import networkx as nx
import ast

class Cohort:
    """Represents a cohort of samples, each with associated variant data."""

    def __init__(self, cohort_df: pd.DataFrame, cohort_df_unfiltered: pd.DataFrame, samples: list, ref: str, workdir: str, family_dict: dict=None):
        
        self.samples = samples
        self.ref = ref
        self.workdir = workdir
        self.cohort_df = cohort_df
        self.cohort_df_unfiltered = cohort_df_unfiltered
        self.missing_variants_df = dict()
        self.dicast_thr = 0.4
        self.family_dict = family_dict


    def _parse_cohort_samples(self, cohort_samples_str):
        """Parse cohort_samples string into list of dicts."""
        if pd.isna(cohort_samples_str):
            return []
        
        if isinstance(cohort_samples_str, list):
            return cohort_samples_str
        
        if isinstance(cohort_samples_str, str):
            try:
                # Try to parse as literal (e.g., "[{'id': 'HG002', 'gt': '0/1'}]")
                return ast.literal_eval(cohort_samples_str)
            except (ValueError, SyntaxError):
                # If that fails, return empty list
                return []
        
        return []

    def get_missing_variants(self):

        cohort_df = self.cohort_df.copy()
        cohort_df.drop_duplicates(subset=['id'], inplace=True)

        # For each sample, find variants where it's not in cohort_samples
        for sample in self.samples:
            # Create a boolean mask for rows where sample is not in cohort_samples
            # cohort_samples now contains list of dicts like [{'id': 'HG002', 'gt': '0/1'}]
            mask = ~cohort_df['cohort_samples'].apply(
                lambda x: any(item.get('id') == sample for item in self._parse_cohort_samples(x))
            )
            
            # Additional family filter: if family_dict is provided, require that cohort_samples 
            # includes at least one sample from the same family
            if self.family_dict is not None:
                if sample not in self.family_dict:
                    raise ValueError(f"Sample '{sample}' not found in family_dict. All samples must be present in family_dict when it is provided.")
                
                family_members = self.family_dict[sample]
                print(cohort_df['cohort_samples'].dtype)
                if family_members:  # Only apply filter if there are family members
                    # Create a mask that checks if any family member is in cohort_samples
                    family_mask = cohort_df['cohort_samples'].apply(
                        lambda x: any(item.get('id') in family_members for item in self._parse_cohort_samples(x))
                    )
                    # Combine with existing mask - variant must be missing from sample AND have family member
                    mask = mask & family_mask
                else:
                    # If no family members, return empty mask (no variants should be considered missing)
                    mask = pd.Series([False] * len(cohort_df), index=cohort_df.index)
            
            # Store the filtered DataFrame in missing_ids
            self.missing_variants_df[sample] = cohort_df[mask].copy()


    def save_missing_variants(self):
        """Save the missing variants to a file."""

        for sample in self.missing_variants_df:
            filename = sample + '_' + self.ref + '.SVs.raw.tsv'
            df = self.missing_variants_df[sample].copy()
            df['sample'] = sample
            print(len(df))
            df.to_csv('/'.join([self.workdir, filename]), sep='\t', index=False, na_rep='NA')


    def load_dicast_predictions(self):
        """Load the dicast predictions for the missing variants."""
        
        dicast_predictions = []
        for sample in self.samples:
            filename = sample + '_' + self.ref + '.SVs.dicast.tsv'
            df = pd.read_csv('/'.join([self.workdir, filename]), sep='\t')
            dicast_predictions.append(df)

        self.dicast_predictions = pd.concat(dicast_predictions, ignore_index=True)


    def update_cohort_information(self):
        """Update cohort_ac and cohort_samples based on dicast predictions."""
        
        # Create a map of variant ids to updated cohort information
        variant_cohort_map = {}
        
        # Process each variant in dicast_predictions
        for _, row in self.dicast_predictions.iterrows():
            variant_id = row['id']
            sample = row['sample']
            dicast_qual = row['dicast_qual']
            
            # If variant is above threshold, include in cohort
            if dicast_qual >= self.dicast_thr:
                if variant_id not in variant_cohort_map:
                    # Find the original variant info
                    orig_var = self.cohort_df[self.cohort_df['id'] == variant_id].iloc[0]
                    cohort_samples_list = self._parse_cohort_samples(orig_var['cohort_samples'])
                    variant_cohort_map[variant_id] = {
                        'samples': [item.get('id') for item in cohort_samples_list],
                        'ac': orig_var['cohort_ac'],
                        'sc': orig_var['cohort_sc'],
                        'gts': [item.get('gt') for item in cohort_samples_list]
                    }
                
                # Add sample if not already in the list
                if sample not in variant_cohort_map[variant_id]['samples']:
                    variant_cohort_map[variant_id]['samples'].append(str(sample))
                    variant_cohort_map[variant_id]['sc'] += 1
                    variant_cohort_map[variant_id]['ac'] += 1
                    variant_cohort_map[variant_id]['gts'].append('0/1')
        self.variant_cohort_map = variant_cohort_map
        
        # Update cohort_df with new cohort information
        self.cohort_df['updated'] = False
        for variant_id, cohort_info in variant_cohort_map.items():
            mask = self.cohort_df['id'] == variant_id
            if mask.any():
                # Create new cohort_samples list format
                new_cohort_samples = [{'id': sample, 'gt': gt} for sample, gt in zip(cohort_info['samples'], cohort_info['gts'])]
                # Assign the same list object to each matching row
                for idx in self.cohort_df.index[mask]:
                    self.cohort_df.at[idx, 'cohort_samples'] = new_cohort_samples
                self.cohort_df.loc[mask, 'cohort_ac'] = cohort_info['ac']
                self.cohort_df.loc[mask, 'cohort_sc'] = cohort_info['sc']
                self.cohort_df.loc[mask, 'updated'] = True

    def find_overlapping_variants(self):
        """Find overlapping variants between cohort_df_unfiltered and dicast_predictions.
        Uses clustering to group overlapping variants and creates sample-specific blacklists for variants to exclude.
        """
        
        def reciprocal_overlap(row: pd.Series) -> float:
            """Compute reciprocal overlap between two intervals for a row in a dataframe."""
            A_start, A_end, B_start, B_end = row['start_1'], row['end_1'], row['start_2'], row['end_2']
            overlap_size = max(0, min(A_end, B_end) - max(A_start, B_start))
            size_A = A_end - A_start
            size_B = B_end - B_start
            
            if size_A == 0 or size_B == 0 or overlap_size <= 0:
                return 0
            
            overlap_A = overlap_size / size_A
            overlap_B = overlap_size / size_B
            
            return min(overlap_A, overlap_B)

        def compute_reciprocal_overlap(df1: pd.DataFrame, df2: pd.DataFrame, sv_type: str, overlap_threshold: float) -> pd.DataFrame:
            """Compute reciprocal overlap between two dataframes."""
            # Filter dataframes by sv_type
            df1_filt = df1[df1['sv_type'] == sv_type].copy().reset_index(drop=True)
            df2_filt = df2[df2['sv_type'] == sv_type].copy().reset_index(drop=True)

            if df1_filt.empty or df2_filt.empty:
                return pd.DataFrame()

            # Make sure dtypes are correct
            df1_filt['start'] = df1_filt['start'].astype(int)
            df1_filt['end'] = df1_filt['end'].astype(int)
            df2_filt['start'] = df2_filt['start'].astype(int)
            df2_filt['end'] = df2_filt['end'].astype(int)
            
            closest_intervals = bf.closest(df1_filt, df2_filt, suffixes=('_1','_2'), k=10)
            closest_intervals = closest_intervals.dropna(subset=['id_1', 'id_2']).reset_index(drop=True)
            
            closest_intervals['reciprocal_overlap'] = closest_intervals.apply(reciprocal_overlap, axis=1)
            closest_intervals = closest_intervals[closest_intervals['reciprocal_overlap'] > overlap_threshold].copy().reset_index(drop=True)
            closest_intervals = closest_intervals[closest_intervals['id_1'] != closest_intervals['id_2']].copy().reset_index(drop=True)
            
            return closest_intervals

        def compute_breakpoint_distance(df1: pd.DataFrame, df2: pd.DataFrame, sv_type: str, distance_threshold: int) -> pd.DataFrame:
            """Compute breakpoint distance for INS variants."""
            # Filter dataframes by sv_type
            df1_filt = df1[df1['sv_type'] == sv_type].copy().reset_index(drop=True)
            df2_filt = df2[df2['sv_type'] == sv_type].copy().reset_index(drop=True)

            if df1_filt.empty or df2_filt.empty:
                return pd.DataFrame()

            # Make sure dtypes are correct
            df1_filt['start'] = df1_filt['start'].astype(int)
            df2_filt['start'] = df2_filt['start'].astype(int)
            
            closest_intervals = bf.closest(df1_filt, df2_filt, suffixes=('_1','_2'), k=10)
            closest_intervals = closest_intervals.dropna(subset=['id_1', 'id_2']).reset_index(drop=True)
            
            # Calculate breakpoint distance (distance between start positions)
            closest_intervals['breakpoint_distance'] = abs(closest_intervals['start_1'] - closest_intervals['start_2'])
            closest_intervals = closest_intervals[closest_intervals['breakpoint_distance'] < distance_threshold].copy().reset_index(drop=True)
            closest_intervals = closest_intervals[closest_intervals['id_1'] != closest_intervals['id_2']].copy().reset_index(drop=True)
            
            return closest_intervals

        def compute_sv_clusters(df: pd.DataFrame, df_overlap: pd.DataFrame, id_column: str) -> pd.DataFrame:
            """Compute clusters of SVs based on an overlap dataframe."""
            df = df.copy()
            
            if df_overlap.empty:
                # No overlaps, each variant gets its own cluster
                df['cluster'] = range(len(df))
                return df
            
            # Create a graph from the overlap dataframe
            G = nx.from_pandas_edgelist(df_overlap, 'id_1', 'id_2')
            
            # Find connected components
            components = list(nx.connected_components(G))
            
            # Create a mapping from node to component
            node_to_component = {}
            for i, component in enumerate(components):
                for node in component:
                    node_to_component[node] = i
                    
            # Map the components back to the dataframe
            df['cluster'] = df[id_column].map(node_to_component)
                    
            # For those without a cluster, assign a new unique cluster label
            next_cluster = len(components)
            for idx, row in df.iterrows():
                if pd.isna(row['cluster']):
                    df.at[idx, 'cluster'] = next_cluster
                    next_cluster += 1
                    
            df['cluster'] = df['cluster'].astype(int)
                    
            return df

        # Initialize sample-specific blacklists and statistics
        self.variants_blacklist = {}  # Dictionary with sample keys
        self.old_variants_blacklist = {}  # Track which blacklisted variants are old per sample
        self.new_variants_blacklist = {}  # Track which blacklisted variants are new per sample
        
        # Process each sample separately
        for sample in self.samples:
            
            # Initialize blacklists for this sample
            self.variants_blacklist[sample] = []
            self.old_variants_blacklist[sample] = []
            self.new_variants_blacklist[sample] = []
            
            # Get high-quality dicast predictions for this sample
            sample_dicast_predictions = self.dicast_predictions[
                (self.dicast_predictions['sample'] == sample) &
                (self.dicast_predictions['dicast_qual'] >= self.dicast_thr)
            ].copy()
            
            if sample_dicast_predictions.empty:
                print(f"No high-quality dicast predictions found for sample {sample}.")
                continue

            # Annotate dicast predictions with cohort information from cohort_df
            sample_cohort_df = self.cohort_df[self.cohort_df['id'].isin(sample_dicast_predictions['id'])].copy().reset_index(drop=True)
            sample_cohort_df.drop_duplicates(subset=['id'], inplace=True)
            sample_dicast_predictions = pd.merge(sample_dicast_predictions, sample_cohort_df[['id', 'cohort_ac', 'cohort_samples']], on='id', how='left')
            sample_dicast_predictions = sample_dicast_predictions[['id', 'sample', 'sv_type', 'chrom', 'start', 'end', 'sv_len',
                                                                   'dicast_qual', 'cohort_ac', 'cohort_samples']].copy().reset_index(drop=True)
            
            # Get variants that are relevant for this sample:
            # 1. Variants from cohort_df_unfiltered that contain this sample
            sample_cohort_variants = self.cohort_df_unfiltered[self.cohort_df_unfiltered['sample'] == sample].copy()
            sample_cohort_variants.rename(columns={'qual': 'dicast_qual'}, inplace=True)
            sample_cohort_variants = sample_cohort_variants[['id', 'sample', 'sv_type', 'chrom', 'start', 'end', 'sv_len',
                                                             'dicast_qual', 'cohort_ac', 'cohort_samples']].copy().reset_index(drop=True)
            
            # Store original variant IDs for statistics
            original_variant_ids = set(sample_cohort_variants['id'])
            
            # Combine both dataframes for clustering (only variants relevant to this sample)
            combined_df = pd.concat([sample_cohort_variants, sample_dicast_predictions], ignore_index=True)
            
            if combined_df.empty:
                print(f"No variants to process for sample {sample}.")
                continue
            
            # Process each SV type
            sv_types = ['DEL', 'DUP', 'INV', 'INS']
            all_overlaps = []
            
            for sv_type in sv_types:
                if sv_type == 'INS':
                    # Use breakpoint distance for INS
                    overlaps = compute_breakpoint_distance(
                        combined_df, combined_df, sv_type, 200
                    )
                else:
                    # Use reciprocal overlap for DEL, DUP, INV
                    overlaps = compute_reciprocal_overlap(
                        combined_df, combined_df, sv_type, 0.5
                    )
                
                if not overlaps.empty:
                    all_overlaps.append(overlaps)
            
            if not all_overlaps:
                print(f"No overlapping variants found for sample {sample}.")
                continue
                
            # Combine all overlaps
            combined_overlaps = pd.concat(all_overlaps, ignore_index=True)
            
            # Compute clusters
            clustered_df = compute_sv_clusters(combined_df, combined_overlaps, 'id')
            
            # For each cluster with more than one variant, choose the best one for this sample
            for cluster_id in clustered_df['cluster'].unique():
                cluster_variants = clustered_df[clustered_df['cluster'] == cluster_id]
                
                if len(cluster_variants) <= 1:
                    continue  # Single variant clusters don't need filtering
                
                # Choose variant with best quality for this sample
                best_variant = None
                best_cohort_ac = -1
                best_dicast_qual = -1
                
                for _, variant in cluster_variants.iterrows():

                    cohort_ac = variant.get('cohort_ac', 0)
                    dicast_qual = variant.get('dicast_qual', 0)

                    # Choose based on highest cohort_ac, tie-break with dicast_qual
                    is_better = (cohort_ac > best_cohort_ac or 
                                (cohort_ac == best_cohort_ac and dicast_qual > best_dicast_qual))
                    
                    if is_better or best_variant is None:
                        best_variant = variant
                        best_cohort_ac = cohort_ac
                        best_dicast_qual = dicast_qual

                # Safety check - if we still don't have a best variant, take the first one
                if best_variant is None:
                    print(f"Warning: No best variant found for cluster {cluster_id}, taking first variant")
                    best_variant = cluster_variants.iloc[0]
                
                # Add all variants except the best one to sample blacklist
                for _, variant in cluster_variants.iterrows():
                    if variant['id'] != best_variant['id']:
                        self.variants_blacklist[sample].append(variant['id'])
                        
                        # Track whether it's old or new
                        if variant['id'] in original_variant_ids:
                            self.old_variants_blacklist[sample].append(variant['id'])
                        else:
                            self.new_variants_blacklist[sample].append(variant['id'])


    def update_csv_file(self, csv_file_path):
        """Update a CSV file with new cohort information, separating by sample and updating existing records."""
        
        # Read the CSV file
        try:
            csv_df = pd.read_csv(csv_file_path)
        except Exception as e:
            print(f"Error reading CSV file {csv_file_path}: {e}")
            return

        # Initialize existing variants for each sample
        self.existing_variants = dict()
        for sample in self.samples:
            self.existing_variants[sample] = csv_df[csv_df['SAMPLE'] == sample]['ID'].tolist()
        
        # First pass: update existing variants with new cohort information
        for idx, row in csv_df.iterrows():
            
            # Get sample and variant id
            sample = row['SAMPLE']
            variant_id = row['ID']

            # Skip variants in sample-specific blacklist
            if variant_id in self.variants_blacklist.get(sample, []):
                continue
            
            # Update cohort info if variant is in variant_cohort_map
            if variant_id in self.variant_cohort_map:
                cohort_info = self.variant_cohort_map[variant_id]
                csv_df.at[idx, 'COHORT_AC'] = cohort_info['ac']
                csv_df.at[idx, 'COHORT_SC'] = cohort_info['sc']
                csv_df.at[idx, 'COHORT_SUP_SAMPLES'] = [{'id': sample, 'gt': gt} for sample, gt in zip(cohort_info['samples'], cohort_info['gts'])]

        # Second pass: add missing variants that pass the threshold
        missing_variants = self.dicast_predictions[self.dicast_predictions['dicast_qual'] >= self.dicast_thr].copy().reset_index(drop=True)
        new_rows = []
        for idx, row in missing_variants.iterrows():
            
            # Get sample and variant id
            sample = row['sample']
            variant_id = row['id']
            
            # Skip variants in sample-specific blacklist or already in CSV
            if variant_id in self.variants_blacklist.get(sample, []) or variant_id in self.existing_variants[sample]:
                continue

            orig_row = csv_df[csv_df['ID'] == variant_id].copy().iloc[0]
            orig_row['SAMPLE'] = sample
            orig_row['QUAL'] = np.round(row['dicast_qual'], 3)
            orig_row['FILTER'] = ['PASS']
            orig_row['GT'] = '0/1'
            orig_row['ID'] = orig_row['SAMPLE'] + '.MERGED.' + orig_row['TYPE'] + '.' + orig_row['CHR'].replace('chr', '') + '.' + orig_row['START'].astype(str) + '.' + orig_row['SIZE'].astype(str)
            orig_row['NUM_SUPP_CALLERS'] = 1
            orig_row['DICAST'] = True
            orig_row['DELLY'] = False
            orig_row['MANTA'] = False
            orig_row['LUMPY'] = False
            orig_row['GRIDSS'] = False
            orig_row['CNVNATOR'] = False
            orig_row['SNIFFLES'] = False
            new_rows.append(orig_row)
        
        # Add new rows to CSV
        csv_df = pd.concat([csv_df, pd.DataFrame(new_rows)], ignore_index=True)
        csv_df.sort_values(by=['SAMPLE', 'ID'], inplace=True)
        
        # Filter out current sample from COHORT_SUP_SAMPLES for all rows
        for idx, row in csv_df.iterrows():
            sample = row['SAMPLE']
            # Parse the supporting samples (handle both string and list formats)
            supporting_samples = row['COHORT_SUP_SAMPLES']
            if isinstance(supporting_samples, str):
                try:
                    supporting_samples = self._parse_cohort_samples(supporting_samples)
                except:
                    continue
            
            if isinstance(supporting_samples, list):
                # Filter out the current sample
                filtered_samples = [entry for entry in supporting_samples if entry.get('id') != sample]
                csv_df.at[idx, 'COHORT_SUP_SAMPLES'] = filtered_samples
        
        csv_df.to_csv(csv_file_path[:-4] + '.regenotyped.csv', sep='\t', index=False)