import pandas as pd

# Load MiRNA CSV
input_file = "/homes/hsuleman/PharmacoGx/PharmacoSet_Exports/NCI60_2021/NCI60_2021_molecular_mirna.csv"
df = pd.read_csv(input_file, index_col=0)

# Transpose data
df_t = df.transpose()

# Save to TSV
output_file = "/homes/hsuleman/PharmacoGx/PharmacoSet_Exports/NCI60_2021/NCI60_mirna.tsv"
df_t.to_csv(output_file, sep='\t')

print(f"Converted TSV saved to: {output_file}")
