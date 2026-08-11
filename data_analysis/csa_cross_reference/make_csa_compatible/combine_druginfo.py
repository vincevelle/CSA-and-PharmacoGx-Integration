import pandas as pd

# Load both TSV files
file1 = pd.read_csv('with_pubchem_ids_noSalt.tsv', sep='\t')
file2 = pd.read_csv('drug_info.tsv', sep='\t')

# Standardize column names for matching
file1 = file1.rename(columns={
    'PubChem_ID': 'PUBCHEM_ID'
})

# Ensure consistent types for comparison
file1['PUBCHEM_ID'] = file1['PUBCHEM_ID'].astype(str)
file2['PUBCHEM_ID'] = file2['PUBCHEM_ID'].astype(str)

# Treat 'nan' strings as actual NaN
file1['PUBCHEM_ID'].replace('nan', pd.NA, inplace=True)
file2['PUBCHEM_ID'].replace('nan', pd.NA, inplace=True)

# Create sets of valid (non-NaN) identifiers
valid_pubchem_ids = set(file1['PUBCHEM_ID'].dropna())
valid_cansmiles = set(file1['canSMILES'].dropna())

# Identify duplicates in file2 with VALID identifiers only
duplicates = file2[
    (file2['PUBCHEM_ID'].isin(valid_pubchem_ids)) |
    (file2['canSMILES'].isin(valid_cansmiles))
]

# Filter out duplicates from file2
file2_filtered = file2.drop(duplicates.index)

# Add missing columns to file1 to match file2 format
for col in file2.columns:
    if col not in file1.columns:
        file1[col] = ''

# Ensure columns are in the same order
file1 = file1[file2.columns]

# Combine the two dataframes
combined = pd.concat([file1, file2_filtered], ignore_index=True)

# Save outputs
combined.to_csv('finalized_druginfo.tsv', sep='\t', index=False)
duplicates.to_csv('duplicates_removed.tsv', sep='\t', index=False)

print("Files merged successfully into 'finalized_druginfo.tsv'")
print("Removed duplicates saved to 'duplicates_removed.tsv'")
