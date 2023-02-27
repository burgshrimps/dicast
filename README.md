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

