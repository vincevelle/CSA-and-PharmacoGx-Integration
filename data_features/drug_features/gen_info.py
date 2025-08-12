import pandas as pd
import pubchempy as pcp
from tqdm import tqdm

# Load the TSV file
df = pd.read_csv("/homes/hsuleman/PharmacoGx/analysis/scripts/drug_SMILES_noSalt.tsv", sep="\t")

failed_smiles = []  # Store SMILES with no CID found

def get_pubchem_cid(smiles):
    """Get the PubChem CID for a given SMILES string using PubChemPy."""
    try:
        compounds = pcp.get_compounds(smiles, namespace='smiles')
        if compounds:
            return compounds[0].cid
        else:
            failed_smiles.append(smiles)
    except Exception as e:
        print(f"Error with SMILES {smiles}: {e}")
        failed_smiles.append(smiles)
    return None

# Fetch PubChem IDs
pubchem_ids = []
for smiles in tqdm(df["canSMILES"], desc="Fetching PubChem IDs"):
    cid = get_pubchem_cid(smiles)
    pubchem_ids.append(cid)

# Add the new column
df["PubChem_ID"] = pubchem_ids

# Save the result
df.to_csv("with_pubchem_ids_noSalt.tsv", sep="\t", index=False)

# Save failed SMILES to a separate file
if failed_smiles:
    pd.Series(failed_smiles).to_csv("failed_smiles.tsv", sep="\t", index=False, header=["SMILES"])
    print(f"{len(failed_smiles)} SMILES could not be resolved. Saved to failed_smiles.tsv")
else:
    print("All SMILES resolved successfully!")

