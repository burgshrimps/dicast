python3 dicast.py cohort \
    --cohort 1kg \
    --ref hg38 \
    --technology ill \
    --workdir /confidential/tGenVar/scripts/tGenVar/dicast/test_cohort \
    --fai /confidential/tGenVar/ref/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai \
    --repeats /confidential/tGenVar/ref/hg38/annotation/hg38_repeatmasker.tsv \
    --cgis /confidential/tGenVar/ref/hg38/annotation/hg38_cpg_islands.tsv \
    --centromeres /confidential/tGenVar/ref/hg38/annotation/hg38_centromeres.tsv \
    --gaps /confidential/tGenVar/ref/hg38/annotation/hg38_asmb_gaps.tsv \
    --althaps /confidential/tGenVar/ref/hg38/annotation/hg38_alt_haps.tsv \
    --vntrs /confidential/tGenVar/ref/hg38/annotation/hg38_vntrs_chaisson.bed \
    --strs /confidential/tGenVar/ref/hg38/annotation/hg38_strs_chaisson.bed \
    --gc /confidential/tGenVar/ref/hg38/annotation/hg38_gc_content.bw \
    --models /confidential/FamilyR13/DATA/10x/SCRIPT/tGenVar/dicast/models/deployment \
    --threads 60 \
    --bams /confidential/FamilyR13/DATA/10x/SCRIPT/tGenVar/dicast/test_cohort/HG001/HG001_GRCh38.bwa_mem.pe.sorted.mdup.recal.bam \
           /confidential/FamilyR13/DATA/10x/SCRIPT/tGenVar/dicast/test_cohort/HG002/HG002_GRCh38.bwa_mem.pe.sorted.mdup.recal.bam \
    --filter-fam \
    --csv /confidential/FamilyR13/DATA/10x/SCRIPT/tGenVar/dicast/test_cohort/Annotated_all.merged.formatted.cohort_ac.csv \
    --ped /confidential/FamilyR13/DATA/10x/SCRIPT/tGenVar/dicast/test_cohort/pedigree.ped
        
