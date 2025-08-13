# Analyzing Data from PharmacoGx

Scripts for visualizing data distributions, calculating summary statistics, and correlating with previous datasets to ensure alignment. Assumes PharmacoDB data has been downloaded into a folder of CSVs using the scripts in /data_retrieval

## 📂 Repository Structure

1. **viability_histograms.py** 
- Loops through the directory of PharmacoDB data, extracts raw viabilities, and plots a histogram contained to the range [0, 200]. Any data outside this range is saved to a separate log file. Find the folder of histograms on lambda7 at /homes/hsuleman/PharmacoGx/data_analysis/viability_histograms
2. **avg_viabilities.py** 
- Calculates the average raw viability values for experiments where the AUC was computed to be > 1. The goal is to find rare cases where the computed AUC is high but the treatment may actually be effective.

### Cross-Referencing With CSA Data

These scripts examine the correlation between the newly computed data from PharmacoDB and the previous response data in the CSA dataset (for experiments that overlap between the two datasets). The goal is to sanity check the curve fitting process and ensure alignment with CSA data

1. **response_correlations.py** 
- Computes correlations and error metrics between overlapping experiments of two response datasets
2. **response_distributions.py** 
- Computes summary statistics and plots distribution of response data as a histogram

#### Make CSA Compatible

Since the sample and drug IDs used in the new PharmacoDB data and the old CSA data are different, we must align ID systems to find experiments that exist in both datasets

1. **combine_druginfo.py** 
- Combines drug_info.tsv files from both PharmacoDB and CSA data, deduplicating based off canonical SMILES and PubChem ID while keeping the PharmacoDB IDs and logging duplicates to a separate file
2. **convert_csa_ids** 
- Uses the combined_druginfo file and the sample_info.csv file to convert the CSA response file to use PharmacoDB IDs, allowing for the matching of overlapping experiments. Note that the sample ID matches from this process are not perfect, since it matches to PharmacoDB ID via cell line name from the DepMap sample_info file, which might not necessarily be exactly the same, so some mappings were almost certainly missed. However for the purposes of validating our response data, we were still able to identify over 100,000 overlapping experiments using this method.
3. **sample_info.csv**
- A sample metadata file downloaded from DepMap (https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap+Public+22Q2&filename=sample_info.csv) containing mappings between DepMap IDs (ACH-XXXXXX) and cell line names. 