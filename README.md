# Dicast

Dicast is a machine learning-based ensemble structural variant caller for short-read sequencing data. The model takes as input a list of variants and outputs the probability of a variant being a true positive call. A typical use case is applying a number of different variant callers on sequencing data of a sample and then using dicast to decide which variants to keep and which ones to discard.

Modes
----------------
There are five modes of operation for dicast.
1. **Preparation:** Reading in the variant calls of all input callers and performs feature extraction.
2. **Training:** Training a new model.
3. **Prediction:** Use an existing model to make predictions on new data.
4. **Evaluation:** Test an existing model.
5. **Manual Curation:** Iteratively train and test an existing model to select a subset of variants for manual curation.

Usage
----------------
Dicast can be executed in the following way\
\
`python dicast.py <mode> <arguments>` \
\
It is necessary to run the preparation mode before being able to run dicast in any other mode.

### Preparation
**Command**\
`python dicast.py [-h] prepare sample ref params workdir`\
\
**Arguments**\
Argument | Description | Optional
--- | --- | ---
sample | Sample name | No
ref | Reference genome name (currently only hg38 supported) | No
params | Preparation parameter file (see below) | No
workdir | Working and output directory. Dicast will create a folder called ensemble in this directory and save the results of the preparation step there. | No \
**Parameter File**\
```javascript
"vcf" : 
    {
        "ill" :
        {
            "delly" : "/confidential/tGenVar/tech/illumina/snakemake_results/sv_REF/delly/SAMPLE/variants.bcf",
            "manta" : "/confidential/tGenVar/tech/illumina/snakemake_results/sv_REF/manta/SAMPLE/results/variants/diploidSV_wINV.vcf.gz",
            "lumpy" : "/confidential/tGenVar/tech/illumina/snakemake_results/sv_REF/lumpy/SAMPLE/SAMPLE-smoove.genotyped.vcf.gz"
        }
    },

    "bam" : 
    {
        "ill": "/confidential/tGenVar/tech/illumina/bam_REF/ill.SAMPLE.REF.bam"
    }     
```

### Training
`python dicast.py [-h] [--chr_exclude CHR] [--cur] svtype params clfparams clfname workdir`
Argument | Description | Optional
--- | --- | ---
svtype | SV type to train the model for | No
params | Training parameter file (see below) | No
clfparams | Model parameter file (see below) | No
clfname | Name of the model in the parameter file | No
chr_excl | List of chromosomes to exclude from training separated by space | Yes
cur | Use this flag to change the labels of variants from manual curation | Yes

