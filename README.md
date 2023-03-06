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
### General
Dicast can be executed in the following way\
\
`python dicast.py <mode> <arguments>` \
\
It is necessary to run the preparation mode before being able to run dicast in any other mode. Example parameter files can be found in `params/`.

### Preparation
This step reads in the VCF files specified in the parameter file and saves their content in the file `<sample>_<ref>.SVs.raw.tsv` inside a directory called `ensemble` within the specified working directory. Features in regards to the genomic context of a variant are saved in the file `<sample>_<ref>.SVs.ref.tsv`. Alignment features are collected in parallel and are written to their respective files with the name `<sample>_<ref>.SVs.aln.<tech>.<chr>.tsv`. Lastly, variant, genomic context, and alignment features are combined into a single output file `<sample>_<ref>.SVs.annot.tsv`.\
\
**Command:** `python dicast.py [-h] prepare sample ref params workdir` \
**Example:** `python dicast.py prepare HG002 hg38 params/params_prep.json /path/to/workdir`
Argument | Description | Optional
--- | --- | ---
sample | Sample name | No
ref | Reference genome name (currently only hg38 supported) | No
params | Preparation parameter file (see below) | No
workdir | Working and output directory. Dicast will create a folder called ensemble in this directory and save the results of the preparation step there. | No 

**Parameters:** `params/params_prep.json` \
The JSON entries under *"ref"* describe the locations of the files containing the reference features for the respective reference genome. Under *"vcf"* there exists one entry for every sequencing technology dicast is applied to. Under every sequencing technology entry there exists an SV caller entry determining the location of the corresponding input VCF file. In the same manner, the location of the input BAM files are declared under *"bam"*. In all entries, *SAMPLE* and *REF* will be replaced with the respective values passed as command line arguments.
- Reference Feature Files
  - Reference Name: Name of the reference genome
    - ...
- VCF
  - Technology
    - SV Caller 1: Location of VCF file 
    - SV Caller 2: Location of VCF file
    - ...
- BAM
  - Technology: Location of BAM file


### Training
With this step, a new model is trained for a specific SV type. The model will be saved as PKL file in the following format `<clfname>_<svtype>.pkl` in the directory specified in `clfparams` under the corresponding model name. 

**Command:** `python dicast.py [-h] [--chr_exclude CHR] [--cur] train svtype params clfparams clfname` \
**Example:** `python dicast.py train DEL params/params_train.json params/params_model.json RF_100`
Argument | Description | Optional
--- | --- | ---
svtype | SV type to train the model for | No
params | Training parameter file (see below) | No
clfparams | Model parameter file (see below) | No
clfname | Name of the model in the parameter file | No
chr_excl | List of chromosomes to exclude from training separated by space | Yes
cur | Use this flag to change the labels of variants from manual curation, if they were falsely labelled  | Yes

**Parameters:** `params/params_train.json` 
- Cohort Name: Describes a group of samples whose file struture is the same
  - Samples: List of samples from this cohort to use for training
  - Reference: Name of the reference genome
  - Variant Features: Output of dicast prepare



