# Downloading Data from PharmacoGx

## 📂 Repository Structure
1. **download_data.R** 
- Downloads datasets from PharmacoDB into a folder of CSV files, organized by sub-dataset
2. **summarize_gene_expr.R** 
- Downloaded rnaseq data comes keyed by SRR ID, this script summarizes data into one gene expression profile per cell line, taking the average (mean) of any replicates


## Available Datasets
1. BeatAML
2. CCLE 
3. CTRPv2 
4. FIMM
5. GCSI
6. GDSCv1/GDSCv2
7. NCI60
8. PRISM
9. PDTX

Find the full folder of downloaded datasets on lambda7 at /homes/hsuleman/PharmacoGx/PharmacoSet_Exports
