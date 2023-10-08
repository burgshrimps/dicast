python3 dicast.py call \
    --cohort 1kg \
    --sample HG002 \
    --ref hg38 \
    --technology ill \
    --workdir /confidential/tGenVar/scripts/tGenVar/dicast/tmp \
    --fai /confidential/tGenVar/ref/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai \
    --repeats /confidential/tGenVar/ref/hg38/annotation/hg38_repeatmasker.tsv \
    --cgis /confidential/tGenVar/ref/hg38/annotation/hg38_cpg_islands.tsv \
    --centromeres /confidential/tGenVar/ref/hg38/annotation/hg38_centromeres.tsv \
    --gaps /confidential/tGenVar/ref/hg38/annotation/hg38_asmb_gaps.tsv \
    --althaps /confidential/tGenVar/ref/hg38/annotation/hg38_alt_haps.tsv \
    --vntrs /confidential/tGenVar/ref/hg38/annotation/hg38_vntrs_chaisson.bed \
    --strs /confidential/tGenVar/ref/hg38/annotation/hg38_strs_chaisson.bed \
    --gc /confidential/tGenVar/ref/hg38/annotation/hg38_gc_content.bw \
    --bam /confidential/tGenVar/scripts/tGenVar/svdb/data/raw/1kg/hg38/HG002/ill/bam/GRCh38.bwa_mem.pe.sorted.mdup.recal.bam \
    --vcfs delly=/confidential/tGenVar/scripts/tGenVar/svdb/data/raw/1kg/hg38/HG002/ill/vcf/delly/formatted_variants.vcf.gz \
           manta=/confidential/tGenVar/scripts/tGenVar/svdb/data/raw/1kg/hg38/HG002/ill/vcf/manta/formatted_variants.vcf.gz \
           lumpy=/confidential/tGenVar/scripts/tGenVar/svdb/data/raw/1kg/hg38/HG002/ill/vcf/lumpy/formatted_variants.vcf.gz \
           cnvnator=/confidential/tGenVar/scripts/tGenVar/svdb/data/raw/1kg/hg38/HG002/ill/vcf/cnvnator/formatted_variants.vcf.gz \
           gridss=/confidential/tGenVar/scripts/tGenVar/svdb/data/raw/1kg/hg38/HG002/ill/vcf/gridss/formatted_variants.vcf.gz \
    --models /confidential/FamilyR13/DATA/10x/SCRIPT/tGenVar/dicast/models/deployment \
    --threads 30 