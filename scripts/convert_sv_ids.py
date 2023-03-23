import pandas as pd
import sys

sample = sys.argv[1]
ref = sys.argv[2]

fname_annot = f'/confidential/FamilyR13/DATA/10x/sv_compare/results/{sample}_{ref}/ensemble/{sample}_{ref}.SVs.annot.tsv' 
fname_curated = f'/confidential/FamilyR13/DATA/10x/sv_compare/results/{sample}_{ref}/curation/{sample}_{ref}.SVs.curated_old_ids.tsv'

df_annot = pd.read_csv(fname_annot, sep='\t', low_memory=False)
df_curation = pd.read_csv(fname_curated, sep='\t', low_memory=False)

df_curation_converted = df_curation.merge(df_annot[['sample', 'tech', 'method', 'type', 'chrom', 'start', 'end', 'id']], on=['sample', 'tech', 'method', 'type', 'chrom', 'start', 'end'], how='left')
df_curation_converted['id'] = df_curation_converted['id_y']
df_curation_converted.drop(['id_x', 'id_y'], axis=1, inplace=True)
print('Dropping:', df_curation_converted['id'].isna().sum())
df_curation_converted.dropna(subset=['id'], inplace=True)

df_curation_converted.to_csv(f'/confidential/FamilyR13/DATA/10x/sv_compare/results/{sample}_{ref}/curation/{sample}_{ref}.SVs.curated.tsv', sep='\t', index=False)