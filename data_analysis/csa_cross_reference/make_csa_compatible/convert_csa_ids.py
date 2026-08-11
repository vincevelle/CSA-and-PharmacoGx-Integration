import pandas as pd

# Load the drug response data
drug_response_df = pd.read_csv("response.tsv", sep="\t")

# Load the chem ID mapping file
chem_map_df = pd.read_csv("duplicates_removed.tsv", sep="\t")  # columns: improve_chem_id, canSMILES, PubChem_ID

# Load the drug info file
finalized_drugs_df = pd.read_csv("finalized_druginfo.tsv", sep="\t")  # columns include: improve_chem_id, canSMILES, PubChem_ID

# Load the sample info file
sample_info_df = pd.read_csv("sample_info.csv")  # columns: DepMap_ID, sample_name

# Create lookup for finalized drug info
canSMILES_to_final = finalized_drugs_df.set_index("canSMILES")["improve_chem_id"].to_dict()
pubchem_to_final = finalized_drugs_df.set_index("PubChem_ID")["improve_chem_id"].to_dict()

# Create lookup for sample info
depmap_to_name = sample_info_df.set_index("DepMap_ID")["cell_line_name"].to_dict()

# Merge chem mapping into drug response to get canSMILES and PubChem_ID for each improve_chem_id
drug_response_df = drug_response_df.merge(
    chem_map_df[["improve_chem_id", "canSMILES", "PubChem_ID"]],
    on="improve_chem_id",
    how="left"
)

# Replace improve_chem_id using finalized_druginfo.tsv
def get_final_chem_id(row):
    if pd.notna(row['canSMILES']) and row['canSMILES'] in canSMILES_to_final:
        return canSMILES_to_final[row['canSMILES']]
    elif pd.notna(row['PubChem_ID']) and row['PubChem_ID'] in pubchem_to_final:
        return pubchem_to_final[row['PubChem_ID']]
    else:
        return row['improve_chem_id']  # No match, keep original

drug_response_df["improve_chem_id"] = drug_response_df.apply(get_final_chem_id, axis=1)

# Replace improve_sample_id using sample_info.csv
drug_response_df["improve_sample_id"] = drug_response_df["improve_sample_id"].map(depmap_to_name).fillna(drug_response_df["improve_sample_id"])

# Drop intermediate columns for clean result
drug_response_df = drug_response_df.drop(columns=["canSMILES", "PubChem_ID"])

# Save result
drug_response_df.to_csv("labeled_response.tsv", sep="\t", index=False)
