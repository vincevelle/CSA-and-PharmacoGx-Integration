import pandas as pd
from rdkit import Chem
from mordred import Calculator, descriptors
from tqdm import tqdm
import multiprocessing as mp
import csv

# Load canonical SMILES
smiles_df = pd.read_csv("/homes/hsuleman/PharmacoGx/analysis/scripts/drug_SMILES_noSalt.tsv", sep="\t")
smiles_df = smiles_df.dropna(subset=["canSMILES"])

# Prepare Mordred calculator
calc = Calculator(descriptors, ignore_3D=True)

# Function to compute descriptors safely
def compute_descriptors(row):
    drug_id = row["improve_chem_id"]
    smiles = row["canSMILES"]
    error_log = []

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")

        descs = calc(mol)
        desc_dict = {}
        for key, value in descs.items():
            col_name = f"mordred.{key}"
            try:
                # Convert value to float if it's a number, else NaN
                if isinstance(value, (int, float)):
                    desc_dict[col_name] = value
                elif value is None:
                    desc_dict[col_name] = float('nan')
                else:
                    raise ValueError(str(value))
            except Exception as e:
                desc_dict[col_name] = float('nan')
                error_log.append({
                    "improve_chem_id": drug_id,
                    "descriptor": col_name,
                    "error": str(e)
                })

        desc_dict["improve_chem_id"] = drug_id
        return (desc_dict, error_log)

    except Exception as e:
        # Molecule couldn't be parsed at all
        return (None, [{
            "improve_chem_id": drug_id,
            "descriptor": "N/A",
            "error": str(e)
        }])

# Use multiprocessing
with mp.Pool(60) as pool:
    results = list(tqdm(
        pool.imap(compute_descriptors, [row for _, row in smiles_df.iterrows()]),
        total=len(smiles_df),
        desc="Calculating Mordred descriptors"
    ))

# Separate valid results and errors
all_descs = []
all_errors = []

for desc_dict, errors in results:
    if desc_dict:
        all_descs.append(desc_dict)
    all_errors.extend(errors)

# Create DataFrame and save results
mordred_df = pd.DataFrame(all_descs).set_index("improve_chem_id")
mordred_df.columns = mordred_df.columns.str.replace(r'[\n\r\t]', '_', regex=True)

# Optionally zero-fill NaNs
mordred_df = mordred_df.fillna(0)

mordred_df.to_csv("drug_mordred.tsv", sep="\t", quoting=csv.QUOTE_MINIMAL)

# Save error log
if all_errors:
    pd.DataFrame(all_errors).to_csv("mordred_descriptor_errors.tsv", sep="\t", index=False)

print(f"Saved {len(mordred_df)} descriptor rows to 'drug_mordred.tsv'")
if all_errors:
    print(f"Logged {len(all_errors)} descriptor-level errors to 'mordred_descriptor_errors.tsv'")
