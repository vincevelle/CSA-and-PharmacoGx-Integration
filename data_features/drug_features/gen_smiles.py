import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from tqdm import tqdm

INPUT_DIR = "PharmacoGx/PharmacoSet_Exports"
OUTPUT_FILE = "drug_SMILES.tsv"

# Function to clean raw SMILES strings
def clean_smiles(s):
    if pd.isna(s):
        return None
    s = str(s).strip().strip('"').strip("'")
    if s.lower() in ["", "na", "null"]:
        return None
    if s.endswith(','):
        s = s[:-1]
    return s

unique_drugs = {}
invalids = []

print(f"Scanning directory: {INPUT_DIR}")
for folder in os.listdir(INPUT_DIR):
    path = os.path.join(INPUT_DIR, folder)
    if not os.path.isdir(path):
        continue

    druginfo_path = os.path.join(path, f"{folder}_druginfo.csv")
    if not os.path.exists(druginfo_path):
        print(f"Skipping {folder} (missing {folder}_druginfo.csv)")
        continue

    try:
        df = pd.read_csv(druginfo_path, index_col=0)
    except Exception as e:
        print(f"Error reading {druginfo_path}: {e}")
        continue

    if 'smiles' not in df.columns:
        print(f"Skipping {folder} (no 'smiles' column)")
        continue

    for drug_id, row in df.iterrows():
        if drug_id not in unique_drugs:
            cleaned = clean_smiles(row['smiles'])
            unique_drugs[drug_id] = cleaned

# Initialize RDKit salt remover
remover = SaltRemover()

# Process and canonicalize SMILES
canonical_smiles = {}
for drug_id, smiles in tqdm(unique_drugs.items(), desc="Canonicalizing SMILES"):
    if smiles is None:
        invalids.append((drug_id, smiles, "null or NA"))
        continue
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            invalids.append((drug_id, smiles, "RDKit parse failed"))
            continue

        # Remove known salts
        mol = remover.StripMol(mol, dontRemoveEverything=True)
        if mol is None or mol.GetNumAtoms() == 0:
            invalids.append((drug_id, smiles, "All fragments removed after salt stripping"))
            continue

        canonical = Chem.MolToSmiles(mol, canonical=True)
        canonical_smiles[drug_id] = canonical
    except Exception as e:
        invalids.append((drug_id, smiles, f"RDKit exception: {e}"))

# Save canonical SMILES
df_out = pd.DataFrame(canonical_smiles.items(), columns=["improve_chem_id", "canSMILES"])
df_out.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"\nSaved {len(df_out)} canonical SMILES to '{OUTPUT_FILE}'")

# Save invalid SMILES for closer inspection
if invalids:
    pd.DataFrame(invalids, columns=["improve_chem_id", "original_smiles", "reason"]) \
        .to_csv("invalid_smiles_canonical.tsv", sep="\t", index=False)
    print(f"Logged {len(invalids)} invalid SMILES to 'invalid_smiles_canonical.tsv'")
