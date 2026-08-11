# Generating Features From PharmacoGx Data

This folder contains scripts for generating features and response metrics from PharmacoDB data. These scripts assume PharmacoDB datasets have been downloaded as CSVs using the scripts from /data_retrieval

## 📂 Repository Structure

### Drug Features

1. **gen_smiles.py** 
- Scans directory of PharmacoDB data to extract a list of valid SMILES strings for each drug. Uses RDKit SaltRemover to strip salts, and then generates canonicalized SMILES. Saves any invalid SMILES to a separate file for closer inspection
2. **gen_mordred.py** 
- Uses the drug_SMILES.tsv file generated from gen_smiles.py to compute Mordred descriptors. Zero fills NaNs and logs errors to a separate file for closer inspection 
1. **gen_fingerprint.py** 
- Uses the drug_SMILES.tsv file generated from gen_smiles.py to compute 512-bit ecfp4 (Extended Connectivity Fingerprint) for all drugs.
2. **gen_info.py** 
- Uses PubChemPy to retrieve PubChem IDs for all canonical SMILES from gen_smiles.py, which can later be used for deduplication when merging with CSA data

### Dose Response Curves
1. **drug_response_curve.py** 
- Called by curve_fitting.py, containing necessary functions for fitting of 4-parameter logistic function to raw experiment data via scipy.optimize.curve_fit(). Uses param_bound = ([low_Einf, low_EC50, low_HS, low_E0], [up_Einf, up_EC50, up_HS, up_E0]) where
param_bound = ([0, -19, -5, 1], [1, 5, 5, 1 + 1e-6]) is initially used, but may be conditionally relaxed to param_bound = ([0, -19, 0, 0], [inf, 5, 5, inf]) in the case where the initial bounds yield a negative HS. See comments in fit_dose_response_curve() function for more detail.
2. **curve_fitting.py** 
- Loops through folder with PharmacoDB data and extracts raw doses and viabilities, which is used to call fit_dose_response_curve() function from drug_response_curve.py to generate response metrics. Multiprocesses curve fitting across 65 processes, adjust based on your machine. 
- Note that raw viabilities are truncated to the range (0,200) prior to curve fitting, this is due to the presence of extreme outliers in some of the datasets (NCI60, for example, has some viabilities over 500). The doses are initially in micromolar (see example plots on https://pharmacodb.ca/), so they are converted to molar for compatibility with drug_response_curve.py




