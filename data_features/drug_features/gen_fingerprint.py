import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.DataStructs import ConvertToNumpyArray
import numpy as np

# Load the file with canonical SMILES
df = pd.read_csv("/homes/hsuleman/PharmacoGx/analysis/scripts/drug_SMILES_noSalt.tsv", sep="\t")

# Create a Morgan fingerprint generator (radius=2 = ECFP4, 512 bits)
generator = GetMorganGenerator(radius=2, fpSize=512)

# Function to compute fingerprint
def smiles_to_ecfp4(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(512, dtype=int)
    fp = generator.GetFingerprint(mol)
    arr = np.zeros((512,), dtype=int)
    ConvertToNumpyArray(fp, arr)
    return arr

# Generate fingerprints
fingerprints = df["canSMILES"].apply(smiles_to_ecfp4)
fp_matrix = np.vstack(fingerprints.values)

# Create the DataFrame
fp_df = pd.DataFrame(fp_matrix,
                     index=df["improve_chem_id"],
                     columns=[f"fp_{i}" for i in range(512)])
fp_df.index.name = ""

# Save to TSV
fp_df.to_csv("drug_fingerprints.tsv", sep="\t")
