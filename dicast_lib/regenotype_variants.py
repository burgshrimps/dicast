import pandas as pd 
import pysam 
import os
import re 
import numpy as np
import logging

class VariantRegenotyping:
    #This class uses the existing class VariantPrep to read all VCF files for one sample. It does this for all samples in the cohort.  
    """ Class to prepare raw variant calls for feature extraction. """

    def __init__(self, cohort: str, ref: str, workdir: str, technology: str, vcfs: list, bams:list , chroms: list, chrom_sizes: str, sv_types: list):
        """ Constructor for VariantPrepCohort class.
        Args:
            cohort (str): Cohort name
            ref (str): Reference genome name
            workdir (str): Working and output directory
            technology (str): Sequencing technology name
            vcfs (list): List of VCF files
            bams (list): List of BAM files
            chroms (list): Chromosomes to use
            chrom_sizes (str): FAI file containing chromosome sizes
            sv_types (list): SV types supported by dicast
        """        

        # Input parameters
        self.cohort = cohort
        self.ref = ref
        self.technology = technology
        self.workdir = workdir
        self.samples_files = self.create_sample_files_dict(vcfs, bams)
        self.samples = self.samples_files.keys()
        self.sv_types = sv_types

        # Auxiallary files for preparation
        self.chrom_sizes = pd.read_csv(chrom_sizes, sep='\t', header=None, 
                                       names=['size', 'offset', 'linebases', 'linewidth'], index_col=0)

        # List of chromosomes to use
        self.chroms = chroms
        
        self.dicast_threshold = 0.4

    def create_sample_files_dict(self, vcfs: list, bams: list) -> dict:
        """ creates a dictionary of samples and their corresponding VCF and BAM files.
        """        
        # get the sample name from the last part of the file path example: calls/merged/GRCh38/HG002/variants.merged.cohort_ac.vcf.gz 
        vcfs_dict = {vcf_path.split('/')[-2]: vcf_path for vcf_path in vcfs}
        bams_dict = {bam_path.split('/')[-2]: bam_path for bam_path in bams}

        # check if both length of the dictionaries are the same
        if len(vcfs_dict.keys()) != len(bams_dict.keys()):
            missing_vcf_files = set(vcfs_dict.keys()) - set(bams_dict.keys())
            missing_bam_files = set(bams_dict.keys()) - set(vcfs_dict.keys())
            if missing_vcf_files:
                raise ValueError(f'VCF files missing for samples: {missing_vcf_files}')
            if missing_bam_files:
                raise ValueError(f'BAM files missing for samples: {missing_bam_files}')
            else:
                raise ValueError('VCF and BAM files do not match.')        
        # mix the two dictionaries by the sample name {sample: (vcf, bam)}
        vcfs_dict = {sample: {'vcf': vcfs_dict[sample], 'bam': bams_dict[sample]} for sample in vcfs_dict.keys()}
        return vcfs_dict
    

    def read_variants_cohort_to_dict(self) :
        """
        Read the records of all vcf files to a dictionary.
        eg. self.records_samples = {'chr1_100_200_DEL': pysam.VariantRecord, 'chr1_300_400_DEL': pysam.VariantRecord}
        """     

        self.records_samples = {}
        self.vcfs_header = pysam.VariantFile(list(self.samples_files.values())[0]['vcf']).header
        for sample in self.samples_files.keys():

            vcf_file = self.samples_files[sample]['vcf']

            if os.path.exists(vcf_file):
                vcf = pysam.VariantFile(vcf_file, "r")
            else:
                # throw error
                raise ValueError('VCF file does not exist.')
            
            
            for rec in vcf.fetch():
                self.records_samples[f"{rec.chrom}_{rec.pos}_{rec.stop}_{rec.info['SVTYPE']}"] = rec
        
            vcf.close()

    def get_missing_varaints_based_on_cohort(self):
        """ get the missing variants based on the cohort.
        """        
        missing_variants_dict = {'id':[], 'ORIGIN_ID':[], 'sv_type':[], 'CALLER':[], 'NUM_SUPP_CALLERS':[], 'end':[], 'start':[],'chrom':[],'chrom_2':[], 'sv_len':[], 'SV_SUBTYPE':[], 'CALLER_Q':[], 'qual':[], 'SUPP_SAMPLES':[], 'SUPP_SAMPLES_GT':[], 'COHORT_AC':[], 'GT':[], 'sv_len':[], 'sample':[], 'cohort':[], 'technology':[], 'caller':[], 'reference':[]} # 'filter':[]
        
        for region in self.records_samples.keys():
            record = self.records_samples[region]
            sup_samples = record.info['SUPP_SAMPLES']

            # get the samples that are not in the supp_samples
            missing_samples = list(set(self.samples) - set(sup_samples))
            current_sample = record.info['SAMPLE_ID']
            for mis_s  in missing_samples:
                missing_variants_dict['ORIGIN_ID'].append(record.info['ORIGIN_ID'])
                missing_variants_dict['sv_type'].append(record.info['SVTYPE'])
                if record.info['SVTYPE'] == 'BND':
                    try:
                        missing_variants_dict['chrom_2'].append(re.search(r'chr.*:', record.alts[0]).group(0)[:-1])
                        missing_variants_dict['end'].append(re.search(r':[0-9]*', record.alts[0]).group(0)[1:])
                    except AttributeError:
                        missing_variants_dict['chrom_2'].append(record.chrom)
                        if record.stop == record.start + 1:
                            missing_variants_dict['end'].append(record.start + 2)
                        else:
                            missing_variants_dict['end'].append(record.stop)
                    missing_variants_dict['sv_len'].append(np.nan)
                else:
                    missing_variants_dict['chrom_2'].append(record.chrom)
                    missing_variants_dict['end'].append(record.stop)
                    missing_variants_dict['sv_len'].append(int(record.info['SVLEN']))
                missing_variants_dict['sample'].append(mis_s)
                missing_variants_dict['cohort'].append('Ashkanazi')
                missing_variants_dict['technology'].append('ill')
                missing_variants_dict['caller'].append('DICAST')
                missing_variants_dict['reference'].append('GRCh38')
                missing_variants_dict['CALLER'].append(record.info['CALLER'])
                missing_variants_dict['NUM_SUPP_CALLERS'].append(record.info['NUM_SUPP_CALLERS'])
                missing_variants_dict['start'].append(record.start)
                missing_variants_dict['chrom'].append(record.chrom)
                missing_variants_dict['SV_SUBTYPE'].append(record.info['SV_SUBTYPE'])
                missing_variants_dict['CALLER_Q'].append(record.info['CALLER_Q'])
                missing_variants_dict['qual'].append(record.qual)
                missing_variants_dict['SUPP_SAMPLES'].append(record.info['SUPP_SAMPLES'])
                missing_variants_dict['SUPP_SAMPLES_GT'].append(record.info['SUPP_SAMPLES_GT'])
                missing_variants_dict['COHORT_AC'].append(record.info['COHORT_AC'])
                # check if the variant is phased 
                gt1 = str(record.samples[current_sample]['GT'][0])
                gt2 = str(record.samples[current_sample]['GT'][1])
                if gt1 == None:
                    gt1 = '.'
                if gt2 == None:
                    gt2 = '.'
                if record.samples[current_sample].phased:
                    missing_variants_dict['GT'].append(str(gt1) + '|' + str(gt2))
                else:
                    missing_variants_dict['GT'].append(str(record.samples[current_sample]['GT'][0]) + '/' + str(record.samples[current_sample]['GT'][1]))

                missing_variants_dict['id'].append(f"{record.chrom}_{record.pos}_{record.stop}_{record.info['SVTYPE']}")

        self.df_variants = pd.DataFrame(missing_variants_dict)
        self.df_variants = self.df_variants[self.df_variants['qual'] > 0.2].reset_index(drop=True)

    def add_dicast_cols(self,missing_variants_features):
        '''Add columns to missing_variants_features dataframe'''
        
        missing_variants_features['merged_id'] =  "NA"
        missing_variants_features['caller_id'] =  "DICAST"
        missing_variants_features['filter'] =  "NA"
        missing_variants_features['caller_qual'] =  "NA"
        missing_variants_features['genotype'] =  missing_variants_features['GT']
        missing_variants_features['performed_confirmation'] =  "NA"
        missing_variants_features['confirmation_status'] =  "NA"
        missing_variants_features['performed_curation'] =  "NA"
        missing_variants_features['curation_status'] =  "NA"
        return missing_variants_features


    def read_regenotyped_variants(self):
        """ Saves variants dataframe to file. """
        all_samples_dfs = []
        for sample in self.df_variants['sample'].unique():
            df = pd.read_csv(f"{self.workdir}/{sample}_{self.ref}.SVs.dicast.tsv" , sep='\t')
            all_samples_dfs.append(df)

        # concatenate the predictions of all samples
        variants_with_dicast_score = pd.concat(all_samples_dfs, ignore_index=True)
   
        # merge based on id start and end 
        self.df_variants_regenotyped = self.df_variants.merge(variants_with_dicast_score, how='left',  on = ['id','cohort', 'sample', 'reference', 'technology', 'caller','sv_type', 'chrom', 'chrom_2', 'start', 'end','sv_len','qual']).reset_index(drop=True)    
        
        # get the ones that are bigger than 0.4
        self.df_variants_regenotyped = self.df_variants_regenotyped[self.df_variants_regenotyped['dicast_qual'] > self.dicast_threshold].reset_index(drop=True)

        # convert the dicast_qual to string
        self.df_variants_regenotyped['dicast_qual'] = self.df_variants_regenotyped['dicast_qual'].astype(str)


        # group the df_variants to get the samples that don't have the variant in one row
        self.df_variants_regenotyped  = self.df_variants_regenotyped.groupby(['id', 'sv_len',  'ORIGIN_ID', 'NUM_SUPP_CALLERS', 'SV_SUBTYPE' , 'SUPP_SAMPLES', 'SUPP_SAMPLES_GT', 'COHORT_AC','CALLER_Q','sv_type']).agg({'sample': ','.join, 'dicast_qual': '|'.join}).reset_index()

        # rename the sample column to SUPP_SAMPLES_REGENOTYPED

        self.df_variants_regenotyped.rename(columns={'sample': 'SUPP_SAMPLES_REGENOTYPED'}, inplace=True)


        # most occuring genotype
        self.df_variants_regenotyped['genotype'] = self.df_variants_regenotyped['SUPP_SAMPLES_GT'].apply(lambda x: max(set(x), key = x.count))

        logging.info(f"# Number of regenotyped variants: {self.df_variants_regenotyped.shape[0]}")
        logging.info(f"# Number of regenotyped DEL variants: {self.df_variants_regenotyped[self.df_variants_regenotyped['sv_type'] == 'DEL'].shape[0]}")
        logging.info(f"# Number of regenotyped INS variants: {self.df_variants_regenotyped[self.df_variants_regenotyped['sv_type'] == 'INS'].shape[0]}")
        logging.info(f"# Number of regenotyped DUP variants: {self.df_variants_regenotyped[self.df_variants_regenotyped['sv_type'] == 'DUP'].shape[0]}")
        logging.info(f"# Number of regenotyped INV variants: {self.df_variants_regenotyped[self.df_variants_regenotyped['sv_type'] == 'INV'].shape[0]}")



        
    def row_to_record(self, row: dict, new_record: pysam.VariantRecord, new_samples_to_add: list, n_new_samples: int, sv_type: str, chrom: str, start: int, end: int, dicast_qual: float, sample: str, additional_gt: str) -> pysam.VariantRecord:
        """
        Converts a row to a VCF record.
        Parameters:
        row (dict): A dictionary containing the row data.
        new_record (vcf.model._Record): A new VCF record object to be populated.
        new_samples_to_add (list): A list of new samples to add.
        n_new_samples (int): The number of new samples.
        sv_type (str): The structural variant type.
        chrom (str): The chromosome.
        start (int): The start position.
        end (int): The end position.
        dicast_qual (float): The quality score.
        sample (str): The sample ID.
        additional_gt (str): Additional genotype information.
        Returns:
        vcf.model._Record: The populated VCF record.
        """
        """converts a row to a vcf record"""

        new_record.chrom = chrom
        new_record.pos = int(start)
        new_record.stop = int(end)
        new_record.id = '.'
        new_record.ref = 'N'
        new_record.alts = [f"<{sv_type}>"]
        new_record.qual = float(dicast_qual)
        new_record.filter.add('PASS')
        
        new_record.info['ORIGIN_ID'] = row['ORIGIN_ID']
        new_record.info['SVTYPE'] = sv_type
        new_record.info['SAMPLE_ID'] = sample
        new_record.info['CALLER'] = 'DICAST'
        new_record.info['NUM_SUPP_CALLERS'] = int(row['NUM_SUPP_CALLERS']) + 1
        new_record.info['SVLEN'] = row['sv_len']
        new_record.info['SV_SUBTYPE'] = row['SV_SUBTYPE']
        new_record.info['CALLER_Q'] = dicast_qual
        new_record.info['SUPP_SAMPLES'] = f"{','.join(list(row['SUPP_SAMPLES']))},{','.join(new_samples_to_add)}"
        new_record.info['SUPP_SAMPLES_GT'] = f"{','.join(list(row['SUPP_SAMPLES_GT']))},{additional_gt}"
        new_record.info['COHORT_AC'] = int(row['COHORT_AC']) + n_new_samples
        for gt in row['SUPP_SAMPLES_GT']:
            if '/' in gt:
                gt_values = gt.split('/')
            elif '|' in gt:
                gt_values = gt.split('|')
   
        new_record.samples[0]['GT'] = (
            int(gt_values[0]) if gt_values[0] != '.' else None, 
            int(gt_values[1]) if gt_values[1] != '.' else None
        )
        
        return new_record

    def load_records_to_dict(self):
        """
        load the vcf records  of the samples to a dictionary
        eg. self.vcf_records = {'HG002': {'chr1_100_200_DEL': pysam.VariantRecord, 'chr1_300_400_DEL': pysam.VariantRecord}, 
                                'HG003': {'chr1_100_200_DEL': pysam.VariantRecord, 'chr1_300_400_DEL': pysam.VariantRecord}}

        """
        self.vcf_records = {}
        for sample in self.samples:
            self.vcf_records[sample] = {}
            vcf_file = self.samples_files[sample]['vcf']
            vcf_in = pysam.VariantFile(vcf_file, "r")
            for rec in vcf_in.fetch():
                self.vcf_records[sample][f"{rec.chrom}_{rec.pos}_{rec.stop}_{rec.info['SVTYPE']}"] = rec
            vcf_in.close()

    def update_variants_and_correct_ac(self):
        """
        Update the variants and correct the allele count (AC) in the VCF records.
        The function iterates over the regenotyped variants and updates the 'SUPP_SAMPLES' 
        and 'COHORT_AC' fields in the VCF records. It also adds new variants to the VCF 
        records if they are not already present.
        The method performs the following steps:
        1. Iterates over each row in the dataframe `self.df_variants_regenotyped`.
        2. Extracts the new samples and their corresponding dicast quality scores.
        3. Updates the existing variant records with new genotyping information.
        4. Adds new variant records to the VCF file if they are not already present..
        """


        # iterate over the regenotyped variants and update the supp_samples and COHORT_AC
        for _, row in self.df_variants_regenotyped.iterrows():
            # samples that have the new variant (these variants have a dicast_qual > 0.4) eg. HG002|HG003
            new_samples_to_add = row['SUPP_SAMPLES_REGENOTYPED'].split(',')

            # corresponding dicast_qual for the new variants for each sample ex. 0.5|0.6 
            dicast_quals_to_add = row['dicast_qual'].split('|')

            # number of new samples
            n_new_samples = len(new_samples_to_add)

            # get the variant id
            var_id = row["id"]
            
            # get the variant info
            chrom, start, end, sv_type = var_id.split('_')

            # gt to add to the supp_samples_gt info field
            additional_gt = ','.join([row['genotype']] * n_new_samples)

            for sample in self.samples:
                # Update existing variant with new genotyping information from other samples
                if var_id in self.vcf_records[sample]: 
                    record = self.vcf_records[sample][var_id]
                    record.info['SUPP_SAMPLES'] = ','.join(list(record.info['SUPP_SAMPLES'])) + "," +  row['SUPP_SAMPLES_REGENOTYPED']
                    record.info['SUPP_SAMPLES_GT'] = ','.join(list(tuple(record.info['SUPP_SAMPLES_GT']))) + "," +  additional_gt
                    record.info['COHORT_AC'] += n_new_samples
                    record.info['NUM_SUPP_CALLERS'] += 1
                # Add new variant to the VCF file 
                elif sample in new_samples_to_add:
                    sample_idx = new_samples_to_add.index(sample)
                    dicast_qual = dicast_quals_to_add[sample_idx]
                    new_record = self.vcfs_header.new_record()
                    rec = self.row_to_record(row, new_record, new_samples_to_add, n_new_samples, sv_type, chrom, start, end, dicast_qual, sample, additional_gt)
                    self.vcf_records[sample][var_id] = rec
                else:
                    continue

    def write_regenotyped_variants(self):
        """
        Write the regenotyped variants to new VCF files after adding new variants and correcting the supp_samples/COHORT_AC.
        This method iterates over the updated records of each sample in `self.samples` and writes the records to new VCF files.
        """

        for sample in self.samples:
            logging.info(f"# Writing regenotyped variants for sample {sample}")
            # old vcf file
            vcf_file = self.samples_files[sample]['vcf']
            
            vcf_file_header = pysam.VariantFile(vcf_file, "r").header

            out_file  = vcf_file.replace('.vcf', '.regenotyped.vcf')
            vcf_out = pysam.VariantFile(out_file, "w", header=vcf_file_header)

            for _, rec in self.vcf_records[sample].items():
                vcf_out.write(rec)
            vcf_out.close()

    def save_missing_variants(self):  
        """ Saves variants dataframe to file. """
        for sample in self.df_variants['sample'].unique():
            self.df_variants[self.df_variants['sample'] == sample].to_csv('/'.join([self.workdir, f"{sample}_{self.ref}.SVs.raw.tsv"]), index=False, sep='\t', na_rep='NA')

