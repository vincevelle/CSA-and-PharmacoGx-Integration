import pandas as pd
import requests
import time
from tqdm import tqdm

# Load the TSV file
df = pd.read_csv("/homes/hsuleman/PharmacoGx/analysis/scripts/drug_SMILES_noSalt.tsv", sep="\t")

def get_pubchem_cid(smiles):
    """Get the PubChem CID for a given SMILES string."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/TXT"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return response.text.strip()
    except Exception as e:
        print(f"Error with SMILES {smiles}: {e}")
    return None

# Fetch PubChem IDs
pubchem_ids = []
for smiles in tqdm(df["canSMILES"], desc="Fetching PubChem IDs"):
    cid = get_pubchem_cid(smiles)
    pubchem_ids.append(cid)

# Add the new column
df["PubChem_ID"] = pubchem_ids

# Save the result to a new TSV file
df.to_csv("with_pubchem_ids_noSalt.tsv", sep="\t", index=False)
