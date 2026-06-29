# CSA Dataset Expansion via PharmacoGx

## 📄 Project Overview
This repository extends the cross-study analysis (CSA) dataset used in the *Drug Response Model Benchmarking* paper (https://arxiv.org/abs/2503.14356v1) by retrieving additional pharmacogenomic data from **PharmacoDB** using the **PharmacoGx** R package. The goal is to align and enrich these datasets to meet the **IMPROVE framework** requirements for running drug response prediction models.

For further information about the processing of multiomics data, see https://github.com/zhuyitan/Data_curation?tab=readme-ov-file 

---

## 📂 Repository Structure
1. **data_retrieval** – Download datasets from PharmacoDB via PharmacoGx.
2. **data_formatting** – Harmonize units, collect metadata, and restructure for IMPROVE compatibility.
3. **data_features** – Compute molecular descriptors, fingerprints, and dose–response metrics.
4. **data_analysis** – Summarize and visualize data distributions and characteristics.

Further details about the scripts for each folder can be found in that folder's README

---

## ⚙️ Installation 

### Environment Setup
```bash
# Create and activate Python environment
conda create -n pgx_env python=3.10
conda activate pgx_env
pip install -r requirements.txt
```

---

## 🔄 Workflow
To reproduce results, start by downloading the data using the scripts provided in ./data_retrieval. This will download PharmacoDB datasets into a directory of CSV files. You may then use the scripts in ./data_features/drug_features to generate molecular features, the scripts in ./data_features/drug_response_curves to fit response curves and compute AUCs, and the scripts in ./data_formatting to convert the raw omics data into CSA units and formatting. 


## References
1. https://github.com/JDACS4C-IMPROVE
2. Smirnov, Petr, et al. "PharmacoDB: an integrative database for mining in vitro anticancer drug screening studies." Nucleic Acids Research (2017).
3. Smirnov, Petr, et al. "PharmacoGx: an R package for analysis of large pharmacogenomic datasets." Bioinformatics 32.8 (2015): 1244-1246.
4. Cross-Study Analysis: https://arxiv.org/abs/2503.14356v1
