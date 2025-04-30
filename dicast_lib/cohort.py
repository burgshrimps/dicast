import pandas as pd
from typing import Dict
import vcfpy
import numpy as np
import os

class Cohort:
    """Represents a cohort of samples, each with associated variant data."""

    def __init__(self, cohort_df: pd.DataFrame, samples: list, ref: str, workdir: str, vcf_files: dict=None):
        
        self.samples = samples
        self.ref = ref
        self.workdir = workdir
        self.cohort_df = cohort_df
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
            self.cohort_df.loc[mask, 'cohort_samples'] = ', '.join(cohort_info['samples'])
            self.cohort_df.loc[mask, 'cohort_ac'] = cohort_info['ac']
            self.cohort_df.loc[mask, 'cohort_samples_gt'] = ', '.join(cohort_info['gts'])
            self.cohort_df.loc[mask, 'updated'] = True

    def update_vcf_files(self):
        
        # Process each sample's VCF
        for sample in self.samples:
            if sample not in self.vcf_files:
                continue
                
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
            
            # First pass: collect existing variants and update their info
            for record in reader:
                variant_id = f"{record.INFO['SVTYPE']}_{record.CHROM}_{record.POS}_{record.INFO['END']}_{record.INFO['SVLEN']}"
 
                existing_variants.add(variant_id)
                # Update cohort info if needed
                if variant_id in self.variant_cohort_map:

                    cohort_info = self.variant_cohort_map[variant_id]
                    record.INFO['COHORT_AC'] = cohort_info['ac']
                    record.INFO['SUPP_SAMPLES'] = cohort_info['samples']
                    record.INFO['SUPP_SAMPLES_GT'] = cohort_info['gts']
                
                updated_records.append(record)
            
            # Write updated existing records
            for record in updated_records:
                writer.write_record(record)
            
            # Second pass: add missing variants that pass the threshold
            missing_variants = self.dicast_predictions[
                (self.dicast_predictions['sample'] == sample) & 
                (self.dicast_predictions['dicast_qual'] >= self.dicast_thr)
            ]
            
            for idx, variant in missing_variants.iterrows():
                variant_id = variant['id']
                
                # Skip if already in VCF
                if variant_id in existing_variants:
                    continue
                    
                # Get full variant information from cohort_df
                var_info = self.cohort_df[self.cohort_df['id'] == variant_id].iloc[0]
                
                # Create a new VCF record
                new_record = self._create_vcf_record(var_info, variant, idx)
                
                if new_record:
                    writer.write_record(new_record)
            
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
                'SUPP_SAMPLES': supp_samples,
                'SUPP_SAMPLES_GT': supp_samples_gt,
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