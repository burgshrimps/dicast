import logging

import pandas as pd

from dicast.vcf_input import (
    VcfInputError,
    canonical_chroms,
    fai_warnings,
    read_caller_vcf,
    read_fai_contigs,
    report_table,
)


class VariantPrep:
    """ Class to prepare raw variant calls for feature extraction. """

    def __init__(self, cohort: str, ref: str, workdir: str, technology: str, chroms: list, chrom_sizes: str, sv_types: list):
        """ Constructor for VariantPrep class.

        Args:
            cohort (str): Cohort name
            ref (str): Reference genome name
            workdir (str): Working and output directory
            chroms (list): Chromosomes to use
            chrom_sizes (str): FAI file containing chromosome sizes
            sv_types (list): SV types supported by dicast
        """

        # Input parameters
        self.cohort = cohort
        self.ref = ref
        self.technology = technology
        self.workdir = workdir
        self.sv_types = sv_types

        # Auxiallary files for preparation
        self.chrom_sizes = pd.read_csv(chrom_sizes, sep='\t', header=None,
                                       names=['size', 'offset', 'linebases', 'linewidth'], index_col=0)

        # List of chromosomes to use, restricted to those the FAI actually
        # has (canonical_chroms raises VcfInputError if none of the
        # requested chromosomes are present at all). fai_warnings flags an
        # FAI that does not look like dicast's chr-named hg38 assumption;
        # neither of these stops the run, they only get logged.
        fai_contigs = read_fai_contigs(chrom_sizes)
        self.chroms, chrom_warnings = canonical_chroms(chroms, fai_contigs)
        for warning in chrom_warnings + fai_warnings(fai_contigs):
            logging.warning(warning)


    def read_vcf(self, vcfs: list, sample: str):
        """ Reads VCF files and stores them in pandas dataframe. """
        self.sample = sample
        self.vcfs = vcfs


    def read_variants(self):
        """ Reads and normalizes every input VCF via
        dicast.vcf_input.read_caller_vcf, one file per (caller, path) pair
        in self.vcfs, and concatenates the results into self.df_variants.

        Raises:
            VcfInputError: if any input file ends up with zero usable
            records (unopenable, unmapped contigs, no matching sample
            column, or every record dropped by the SV contract) -- running
            the rest of the pipeline on such a file would silently drop it
            from every downstream step instead of failing loudly.
        """

        vcf_dfs = []
        self.reports = []
        for caller, vcf_file in self.vcfs:
            df, report = read_caller_vcf(
                vcf_file, caller, self.sample, self.chroms, self.cohort,
                self.ref, self.technology)
            vcf_dfs.append(df)
            self.reports.append(report)

        table = report_table(self.reports)
        for line in table.splitlines():
            logging.info(line)

        zero_usable = [r.path for r in self.reports if r.usable == 0]
        if zero_usable:
            raise VcfInputError(
                'The following input file(s) have zero usable records: '
                f'{", ".join(zero_usable)}\n\n{table}\n\n'
                "Run 'dicast check' on these files for the full report."
            )

        # Merge all VCF files
        self.df_variants = pd.concat(vcf_dfs, ignore_index=True)
        assert self.df_variants['id'].is_unique, 'internal record ids are not unique across input files'


    def check_out_of_bounds(self, svtype: str, chrom: str, chrom_2: str, start: int, end: int, chrom_sizes: pd.DataFrame, padding: int=50) -> bool:
        """ Checks if variant is out of chromosome bounds.

        Args:
            svtype (str): SV type
            chrom (str): Chromosome
            chrom2 (str): Second chromosome for translocations
            start (int): Start position
            end (int): End position
            chrom_sizes (pd.DataFrame): Dataframe with chromosome sizes
            padding (int, optional): Padding around SV borders. Defaults to 50.

        Returns:
            bool: True if variant is out of bounds, False otherwise
        """        

        if svtype != 'BND':
            return start - padding < 0 or end + padding > chrom_sizes.loc[chrom, 'size']
        else:
            # For translocations, check both chromosomes
            outbounds_chrom1 = start - padding < 0 or start + padding > chrom_sizes.loc[chrom, 'size']
            outbounds_chrom2 = end - padding < 0 or end + padding > chrom_sizes.loc[chrom_2, 'size']
            return outbounds_chrom1 or outbounds_chrom2


    def filter_variants(self):
        """ Removes variants that are out of chromosomes bounds or have other problems.

        Raises:
            VcfInputError: if restricting to self.sv_types leaves nothing --
            running the rest of the pipeline on an empty frame would fail
            downstream with a much less informative error.
        """

        # read_variants already rejects runs with nothing usable; this only
        # guards direct callers, since boolean-indexing an empty frame below
        # would fail with an unrelated KeyError.
        if self.df_variants.empty:
            return

        # Remove calls that are out of chromosome bounds
        self.df_variants['start'] = self.df_variants['start'].astype(int)
        self.df_variants['end'] = self.df_variants['end'].astype(int)
        self.df_variants['outbounds'] = self.df_variants.apply(lambda x: self.check_out_of_bounds(x['sv_type'], x['chrom'], x['chrom_2'], x['start'], x['end'], self.chrom_sizes), axis=1)
        self.df_variants = self.df_variants[~self.df_variants['outbounds']].copy().drop('outbounds', axis=1).reset_index(drop=True)

        # Remove SV types that are currently not supported by dicast
        self.df_variants = self.df_variants[self.df_variants['sv_type'].isin(self.sv_types)].copy().reset_index(drop=True)
        if self.df_variants.empty:
            raise VcfInputError(
                f'No variants remain after restricting to SV type(s) {", ".join(self.sv_types)}.'
            )


    def get_variant_df(self):
        return self.df_variants


    def save_variants(self):
        """ Saves variants dataframe to file, under workdir/input/. """

        filename = self.sample + '_' + self.ref + '.SVs.raw.tsv'
        self.df_variants.to_csv('/'.join([self.workdir, 'input', filename]), index=False, sep='\t', na_rep='NA')