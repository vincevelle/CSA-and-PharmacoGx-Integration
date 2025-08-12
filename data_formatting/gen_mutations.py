import os
import pandas as pd

# Path to the gene info mapping file
GENE_INFO_PATH = '/homes/hsuleman/PharmacoGx/analysis/outputs/gene_info.csv'

# Load gene info mapping
gene_info = pd.read_csv(GENE_INFO_PATH, dtype=str)
gene_info = gene_info.dropna(subset=['gene_symbol', 'entrez_id', 'ensembl_id'])
gene_info = gene_info.drop_duplicates(subset='gene_symbol').set_index('gene_symbol')

def process_mutation_file(file_path, output_dir):
    df = pd.read_csv(file_path, index_col=0, dtype=str)
    df.dropna(how='all', axis=0, inplace=True)   # Drop entirely NaN rows
    df.dropna(how='all', axis=1, inplace=True)   # Drop entirely NaN columns

    gene_symbols = df.index.astype(str)
    valid_genes = gene_symbols[gene_symbols.isin(gene_info.index)]
    if valid_genes.empty:
        print(f"Skipping {file_path}: no valid gene mappings found.")
        return

    df = df.loc[valid_genes]

    # Fill NaNs with 0s
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

    # Extract gene metadata
    entrez_ids = gene_info.loc[valid_genes, 'entrez_id'].values
    ensembl_ids = gene_info.loc[valid_genes, 'ensembl_id'].values
    gene_names = valid_genes.values

    # Transpose so samples become rows and genes become columns
    df = df.T
    df.columns = gene_names

    # Create metadata rows
    metadata = pd.DataFrame([entrez_ids, gene_names, ensembl_ids], columns=df.columns)
    metadata.index = ['entrez_id', 'gene_symbol', 'ensembl_id']

    # Add placeholder index for metadata rows
    metadata.insert(0, 'sample_id', [''] * 3)

    # Reset sample index to column for saving
    df.insert(0, 'sample_id', df.index)
    df.reset_index(drop=True, inplace=True)

    # Combine metadata and data
    full_df = pd.concat([metadata, df], ignore_index=True)

    # Save to TSV
    dataset_name = os.path.basename(file_path).replace('_mutation.csv', '')
    output_path = os.path.join(output_dir, f"{dataset_name}_processed.tsv")
    os.makedirs(output_dir, exist_ok=True)
    full_df.to_csv(output_path, sep='\t', index=False, header=False)
    print(f"Processed: {file_path} → {output_path}")

def process_all_mutations(base_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for root, _, files in os.walk(base_dir):
        for file in files:
            # Only process exact matches like CCLE_2015_mutation.csv, use _"mutation_binary" for binary mutation files
            if file.endswith('_mutation_binary.csv'):
                file_path = os.path.join(root, file)
                process_mutation_file(file_path, output_dir)

# Main
if __name__ == "__main__":
    BASE_DIR = '/homes/hsuleman/PharmacoSet_Exports'
    OUTPUT_DIR = 'processed_mutations_tsvs'
    process_all_mutations(BASE_DIR, OUTPUT_DIR)

