import pandas as pd
from typing import Dict
import vcfpy
import numpy as np
import os
import bioframe as bf
import networkx as nx

class Cohort:
    """Represents a cohort of samples, each with associated variant data."""

    def __init__(self, cohort_df: pd.DataFrame, cohort_df_unfiltered: pd.DataFrame, samples: list, ref: str, workdir: str, vcf_files: dict=None):
        
        self.samples = samples
        self.ref = ref
        self.workdir = workdir
        self.cohort_df = cohort_df
        self.cohort_df_unfiltered = cohort_df_unfiltered
        self.vcf_files = vcf_files
        self.missing_variants_df = dict()
        self.dicast_thr = 0.4


    def get_missing_variants(self):

        cohort_df = self.cohort_df.copy()
        cohort_df.drop_duplicates(subset=['id'], inplace=True)

        # For each sample, find variants where it's not in cohort_samples
        for sample in self.samples:
            # Create a boolean mask for rows where sample is not in cohort_samples
            mask = ~cohort_df['cohort_samples'].str.contains(sample, na=False)
            
            # Store the filtered DataFrame in missing_ids
            self.missing_variants_df[sample] = cohort_df[mask].copy()


    def save_missing_variants(self):
        """Save the missing variants to a file."""

        for sample in self.missing_variants_df:
            filename = sample + '_' + self.ref + '.SVs.raw.tsv'
            df = self.missing_variants_df[sample].copy()
            df['sample'] = sample
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
                    variant_cohort_map[variant_id] = {
                        'samples': orig_var['cohort_samples'].split(', ') if isinstance(orig_var['cohort_samples'], str) else [],
                        'ac': orig_var['cohort_ac'],
                        'gts': orig_var['cohort_samples_gt'].split(', ') if isinstance(orig_var['cohort_samples_gt'], str) else []
                    }
                
                # Add sample if not already in the list
                if sample not in variant_cohort_map[variant_id]['samples']:
                    variant_cohort_map[variant_id]['samples'].append(sample)
                    variant_cohort_map[variant_id]['ac'] += 1
                    variant_cohort_map[variant_id]['gts'].append('0/1')
        self.variant_cohort_map = variant_cohort_map
        
        # Update cohort_df with new cohort information
        self.cohort_df['updated'] = False
        for variant_id, cohort_info in variant_cohort_map.items():
            mask = self.cohort_df['id'] == variant_id
            if mask.any():
                self.cohort_df.loc[mask, 'cohort_samples'] = ','.join(cohort_info['samples'])
                self.cohort_df.loc[mask, 'cohort_ac'] = cohort_info['ac']
                self.cohort_df.loc[mask, 'cohort_samples_gt'] = ','.join(cohort_info['gts'])
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
            sample_dicast_predictions = pd.merge(sample_dicast_predictions, sample_cohort_df[['id', 'cohort_ac', 'cohort_samples', 'cohort_samples_gt']], on='id', how='left')
            sample_dicast_predictions = sample_dicast_predictions[['id', 'sample', 'sv_type', 'chrom', 'start', 'end', 'sv_len',
                                                                   'dicast_qual', 'cohort_ac', 'cohort_samples', 'cohort_samples_gt']].copy().reset_index(drop=True)
            
            # Get variants that are relevant for this sample:
            # 1. Variants from cohort_df_unfiltered that contain this sample
            sample_cohort_variants = self.cohort_df_unfiltered[self.cohort_df_unfiltered['sample'] == sample].copy()
            sample_cohort_variants.rename(columns={'qual': 'dicast_qual'}, inplace=True)
            sample_cohort_variants = sample_cohort_variants[['id', 'sample', 'sv_type', 'chrom', 'start', 'end', 'sv_len',
                                                             'dicast_qual', 'cohort_ac', 'cohort_samples', 'cohort_samples_gt']].copy().reset_index(drop=True)
            
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


    def update_vcf_files(self):
        
        statistics_dict = dict()
        
        # Process each sample's VCF
        for sample in self.samples:
            if sample not in self.vcf_files:
                continue

            statistics_dict[sample] = {'variants_total': 0, 
                                       'variants_updated': 0, 
                                       'variants_added': 0,
                                       'variants_dicast': 0,
                                       'variants_old_skipped': 0,
                                       'variants_new_skipped': 0,
                                       'variants_exist_skipped': 0,
                                       'variants_total_ids' : [],
                                       'variants_updated_ids' : [],
                                       'variants_added_ids' : [],
                                       'variants_dicast_ids' : [],
                                       'variants_old_skipped_ids' : [],
                                       'variants_new_skipped_ids' : [],
                                       'variants_exist_skipped_ids' : []}
                
            vcf_path = self.vcf_files[sample]
            vcf_out_path = vcf_path.replace('.vcf', '.regenotyped.vcf')
            
            # Read the VCF
            reader = vcfpy.Reader.from_path(vcf_path)
            
            # Create a writer with the same header
            header = reader.header
            
            writer = vcfpy.Writer.from_path(vcf_out_path, header)
            
            # Create a map of existing variant IDs in this sample's VCF
            existing_variants = set()
            updated_records = []
            
            # Get sample-specific blacklist
            sample_blacklist = self.variants_blacklist.get(sample, [])
            
            # First pass: collect existing variants and update their info
            for record in reader:
                variant_id = f"{record.INFO['SVTYPE']}_{record.CHROM}_{record.POS}_{record.INFO['END']}_{record.INFO['SVLEN']}"
                existing_variants.add(variant_id)
                statistics_dict[sample]['variants_total'] += 1
                statistics_dict[sample]['variants_total_ids'].append(variant_id)
                
                # Skip variants in sample-specific blacklist
                if variant_id in sample_blacklist:
                    statistics_dict[sample]['variants_old_skipped'] += 1
                    statistics_dict[sample]['variants_old_skipped_ids'].append(variant_id)
                    continue
                
                # Update cohort info if needed
                if variant_id in self.variant_cohort_map:
                    cohort_info = self.variant_cohort_map[variant_id]
                    record.INFO['COHORT_AC'] = cohort_info['ac']
                    record.INFO['SUPP_SAMPLES'] = cohort_info['samples']
                    record.INFO['SUPP_SAMPLES_GT'] = cohort_info['gts']
                    statistics_dict[sample]['variants_updated'] += 1
                    statistics_dict[sample]['variants_updated_ids'].append(variant_id)
                
                updated_records.append(record)
            
            # Write updated existing records
            for record in updated_records:
                writer.write_record(record)
            
            # Second pass: add missing variants that pass the threshold
            missing_variants = self.dicast_predictions[
                (self.dicast_predictions['sample'] == sample) & 
                (self.dicast_predictions['dicast_qual'] >= self.dicast_thr)
            ]
            statistics_dict[sample]['variants_dicast'] += len(missing_variants)
            statistics_dict[sample]['variants_dicast_ids'].extend(missing_variants['id'].tolist())
            
            for idx, variant in missing_variants.iterrows():
                variant_id = variant['id']
                
                # Skip if already in VCF
                if variant_id in existing_variants:
                    statistics_dict[sample]['variants_exist_skipped'] += 1
                    statistics_dict[sample]['variants_exist_skipped_ids'].append(variant_id)
                    continue
                
                # Skip variants in sample-specific blacklist
                if variant_id in sample_blacklist:
                    statistics_dict[sample]['variants_new_skipped'] += 1
                    statistics_dict[sample]['variants_new_skipped_ids'].append(variant_id)
                    continue
                    
                # Get full variant information from cohort_df
                var_info = self.cohort_df[self.cohort_df['id'] == variant_id].iloc[0]
                
                # Create a new VCF record
                new_record = self._create_vcf_record(var_info, variant, idx)
                
                if new_record:
                    writer.write_record(new_record)
                    statistics_dict[sample]['variants_added'] += 1
                    statistics_dict[sample]['variants_added_ids'].append(variant_id)

            writer.close()


    def _create_vcf_record(self, var_info, dicast_var, idx):
        """Create a VCF record for a missing variant."""
        
        try:
            # Basic variant fields
            chrom = var_info['chrom']
            pos = int(var_info['start'])
            id = '.'
            ref = 'N'  # Placeholder reference allele
            alt = [vcfpy.SymbolicAllele(var_info['sv_type'])]
            qual = float(dicast_var['dicast_qual']) if 'dicast_qual' in dicast_var else None
            filt = ['PASS']
            
            # Format SUPP_SAMPLES and SUPP_SAMPLES_GT properly
            supp_samples = var_info['cohort_samples'].split(', ') if isinstance(var_info['cohort_samples'], str) else []
            supp_samples_gt = var_info['cohort_samples_gt'].split(', ') if isinstance(var_info['cohort_samples_gt'], str) else []
            
            # INFO fields
            info = {
                'ORIGIN_ID': 'dicast.' + var_info['sv_type'] + '.' + str(idx),
                'SVTYPE': var_info['sv_type'],
                'SAMPLE_ID': dicast_var['sample'],
                'CALLER': 'DICAST',
                'NUM_SUPP_CALLERS': 1,
                'END': int(var_info['end']),
                'SVLEN': int(var_info['sv_len']),
                'SV_SUBTYPE': var_info['sv_type'],
                'CALLER_Q': np.round(float(dicast_var['dicast_qual']), 4),
                'DICAST_Q': np.round(float(dicast_var['dicast_qual']), 4),
                'SUPP_SAMPLES': ','.join(supp_samples),
                'SUPP_SAMPLES_GT': ','.join(supp_samples_gt),
                'COHORT_AC': int(var_info['cohort_ac']),
            }
            
            # Format and samples
            fmt = ['GT']
            calls = [vcfpy.Call(dicast_var['sample'], {'GT': '0/1'})]
            
            return vcfpy.Record(
                CHROM=chrom,
                POS=pos,
                ID=[id],
                REF=ref,
                ALT=alt,
                QUAL=qual,
                FILTER=filt,
                INFO=info,
                FORMAT=fmt,
                calls=calls
            )
            
        except Exception as e:
            print(f"Error creating VCF record for {var_info['id']}: {e}")
            return None