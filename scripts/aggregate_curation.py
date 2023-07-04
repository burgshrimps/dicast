import sys 
import pandas as pd 
from glob import glob


def consensus_curate(votes):
    """ Performs sort of majority voting on the curation votes.
    
    param votes: A pandas series of votes. 1 for a positive vote, 0 for a negative vote, and NaN for no vote. 
    
    return 1 if the majority of votes are positive, 0 otherwise. """
    
    num_votes = votes.notna().sum()
    if votes.sum() >= num_votes / 2:
        return 1
    else:
        return 0


def aggregate_curation(sample, ref, date, directory, curators):
    """ Aggregates the curation votes for each chunk and does consensus voting.

    param sample: The sample name.
    param date: The date of the curation.
    param directory: The directory where the curation votes are stored.
    param curators: The list of curators. """
    
    curation_dfs = []
    curator_cols = []
    for curator in curators:
        curator_col = 'Confirmed (' + curator + ')'
        try:
            curation_df = pd.concat([pd.read_csv(fname, sep='\t') for fname in glob(directory + '/*/*/*/*/*' + curator + '_curated.tsv')], ignore_index=True)
            curation_df = curation_df[['sample', 'tech', 'method', 'id', 'type', 'chrom', 'chrom2', 'start', 'end', 'confirmed', 'pred_dicast', 'err_type', curator_col]].copy()
            curation_df[curator_col] = curation_df[curator_col].replace({True: 1, False: 0})
            curation_dfs.append(curation_df)
            curator_cols.append(curator_col)
        except ValueError:
            pass

    curation_df = curation_dfs[0]
    if len(curation_dfs) > 1:
        for i in range(1, len(curation_dfs)):
            curation_df = curation_df.merge(curation_dfs[i], on=['sample', 'tech', 'method', 'id', 'type', 'chrom', 'chrom2', 'start', 'end', 'confirmed', 'pred_dicast', 'err_type'], how='left')
        curation_df['Confirmed (Consensus)'] = curation_df[curator_cols].apply(consensus_curate, axis=1)
    else:
        curation_df['Confirmed (Consensus)'] = curation_df[curator_cols[0]]
        
    curation_df.to_csv(directory + '/' + sample + '_' + ref + '.SVs.curated.tsv', sep='\t', index=False)


date = sys.argv[1]
ref = 'hg38'
samples = ['17-08', '2048-07', '906-02', '18-5047']
curators = ['Nico']
directory = '/confidential/tGenVar/scripts/tGenVar/dicast/training/curation/' + ref

for sample in samples:
    aggregate_curation(sample, ref, date, directory + '/' + sample, curators)

