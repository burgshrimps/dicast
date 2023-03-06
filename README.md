# Dicast

Dicast is a machine learning-based ensemble structural variant caller for short-read sequencing data. The model takes as input a list of variants and outputs the probability of a variant being a true positive call. A typical use case is applying a number of different variant callers on sequencing data of a sample and then using dicast to decide which variants to keep and which ones to discard.

## Modes
There are five modes of operation for dicast.
1. **[Preparation:](#preparation)** Reading in the variant calls of all input callers and performs feature extraction.
2. **[Training:](#training)** Training a new model.
3. **[Prediction:](#prediction)** Use an existing model to make predictions on new data.
4. **[Evaluation:](#evaluation)** Test an existing model.
5. **[Manual Curation:](#manual-curation)** Iteratively train and test an existing model to select a subset of variants for manual curation.

## General
Dicast can be executed in the following way\
\
`python dicast.py <mode> <arguments>` \
\
It is necessary to run the preparation mode before being able to run dicast in any other mode. Example parameter files can be found in `params/`.

## Preparation
This step reads in the VCF files specified in the parameter file and saves their content in the file `<sample>_<ref>.SVs.raw.tsv` inside a directory called `ensemble` within the specified working directory. Features in regards to the genomic context of a variant are saved in the file `<sample>_<ref>.SVs.ref.tsv`. Alignment features are collected in parallel and are written to their respective files with the name `<sample>_<ref>.SVs.aln.<tech>.<chr>.tsv`. Lastly, variant, genomic context, and alignment features are combined into a single output file `<sample>_<ref>.SVs.annot.tsv`.\
\
**Command:** `python dicast.py prepare [-h] sample ref params workdir` \
**Example:** `python dicast.py prepare HG002 hg38 params/params_prep.json /path/to/workdir`
Argument | Description | Optional
--- | --- | ---
sample | Sample name | No
ref | Reference genome name (currently only hg38 supported) | No
params | Preparation parameter file (see below) | No
workdir | Working and output directory. Dicast will create a folder called ensemble in this directory and save the results of the preparation step there. | No 

**Parameters:** `params/params_prep.json` \
The JSON entries under *"ref"* describe the locations of the files containing the reference features for the respective reference genome. Under *"vcf"* there exists one entry for every sequencing technology dicast is applied to. Under every sequencing technology entry there exists an SV caller entry determining the location of the corresponding input VCF file. In the same manner, the location of the input BAM files are declared under *"bam"*. In all entries, *SAMPLE* and *REF* will be replaced with the respective values passed as command line arguments.

## Training
With this step, a new model is trained for a specific SV type. The model will be saved as PKL file in the following format `<clfname>_<svtype>.pkl` in the directory specified in `clfparams` under the corresponding model name. 

**Command:** `python dicast.py train [-h] [--chr_exclude CHR] [--cur] svtype params clfparams clfname` \
**Example:** `python dicast.py train --cur DEL params/params_train.json params/params_model.json RF_100`
Argument | Description | Optional
--- | --- | ---
svtype | SV type to train the model for | No
params | Training parameter file (see below) | No
clfparams | Model parameter file (see below) | No
clfname | Name of the model in the parameter file | No
chr_excl | List of chromosomes to exclude from training separated by space (default: none) | Yes
cur | Use this flag to change the labels of variants from manual curation, if they were falsely labelled  | Yes

**Parameters:** `params/params_train.json` \
The training parameter file consists of one or multiple cohort entries. The example parameter file contains one cohort called "tgenvar". A cohort describes a group of samples with the same directory and file structure. The entry *"samples"* is a list of samples used for training a model, *"ref"* is the name of the reference genome, *"variant_features"* is the output from the preparation step, *"variant_labels"* contains the IDs of all confirmed (true positive) variants. Lastly *"variant_curation"* is the output file of the manual curation and contains corrected labels for manually curated SVs. In all entries, *SAMPLE* and *REF* will be replaced with the respective values from *"samples"* and *"ref"*. \
\
**Parameters:** `params/params_model.json` \
This JSON file contains entries for different models, specified by their name. A corresponding model is selected by specifying the name in the command line argument `clfname`. The entry *"classifier"* describes the type of classifier to train. The entry *"directory"* describes the directory the model is saved in. The entry *"parameters"* describes a number of model-specific parameters.

## Prediction
This step predicts whether a variant is a true positive or a false positive call for a set of unseen SVs. By default the predictions are saved in a file called `<sample>.dicast.tsv` in the `workdir`. The column "qual_dicast" corresponds to the probability of a variant being a true positive call. Additionally, if specified, an INFO tag, DQ, is added to the input VCF files, which corresponds to the probability of a variant being a true postive call. These modified VCF files are saved with the suffix `.dicast.vcf` in the directories of the input VCF files. It is necessary to have a trained model in the model directory for all SV types one wants to make predictions about. Furthermore,  it is necessary to run dicast prepare before running dicast predict.

**Command:** `python dicast.py predict [-h] [--vcf] params clfparams clfname workdir` \
**Example:** `python dicast.py predict --vcf params/params_predict.json params/params_model.json RF_100 /path/to/workdir`
Argument | Description | Optional
--- | --- | ---
params | Prediction parameter file (see below) | No
clfparams | Model parameter file (see below) | No
clfname | Name of the model in the parameter file | No
workdir | Working and output directory | No 
vcf | Additionally add an INFO tag to the input VCF files and save them in workdir | Yes

**Parameters:** `params/params_predict.json` \
The entry *"samples"* specifies a list of samples one wants to apply dicast to. The entry *"variant_features"* is the output of dicast prepare. Under the entry *"vcf"* the input VCF files are specified in case one wants to add the INFO tag to them. In all entries, *SAMPLE* and *REF* will be replaced with the respective values from *"samples"* and *"ref"*. \
\
**Parameters:** `params/params_model.json` \
See above.

## Evaluation
This step uses trained models to make predicitons for a set of unseen SVs and subsequently compares the predictions with the ground truth labels. This results of this mode are saved in a directory called `eval` inside the model directory. For the specified SV type model a file, `<clfname>_<svtype>_eval.tsv` is created containing the feature values, quality scores of the input callers as well as the dicast probabilities. This file can be used by the script `scripts/evaluate_model.ipynb` to compute metrics such as precision and recall.

**Command:** `python dicast.py test [-h] [--chr_incl] [--cur] svtype params clfparams clfname` \
**Example:** `python dicast.py test --cur DEL params/params_test.json params/params_model.json RF_100`
Argument | Description | Optional
--- | --- | ---
svtype | SV type to train the model for | No
params | Training parameter file (see below) | No
clfparams | Model parameter file (see below) | No
clfname | Name of the model in the parameter file | No
chr_incl | List of chromosomes to include in testing separated by space (default: all) | Yes
cur | Use this flag to change the labels of variants from manual curation, if they were falsely labelled  | Yes

**Parameters:** `params/params_test.json` \
Same format as `params/params_train.json`. \
\
**Parameters:** `params/params_model.json` \
See above.

## Manual Curation
For a given set of samples and an SV type, this step iteratively trains a model on all but one chromosome and tests the model on the remaining chromosome. Next, inconsistencies between the ground truth labels and the dicast predictions are identified. There are two kinds of inconsistencies:
1. False Positives (FP): Ground truth says a variant is not confirmed but dicast says it is confirmed
2. False Negatives (FN): Ground truth says a variant is confirmed but dicast says it is not confirmed



