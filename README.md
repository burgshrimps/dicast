# Dicast

Dicast is a machine learning-based ensemble structural variant caller for short-read sequencing data. 

Modes
----------------
There are five modes of operation for Dicast.
1. Preparation: Reading in the variant calls of all input callers and performs feature extraction.
2. Training: Training a new model.
3. Prediction: Use an existing model to make predictions on new data.
4. Evaluation: Test an existing model.
5. Manual Curation: Iteratively train and test an existing model to select a subset of variants for manual curation.

