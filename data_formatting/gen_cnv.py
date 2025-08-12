import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ====== File paths ======
EXPR_CSV = "/homes/hsuleman/PharmacoSet_Exports3/CCLE_2015/CCLE_2015_molecular_cnv.csv"   # First col = Feature (gene symbol), rest = SNP.arrays IDs
CELLINFO_CSV = "/homes/hsuleman/PharmacoSet_Exports3/CCLE_2015/CCLE_2015_cellinfo.csv"     # Contains mappings between 'SNP.arrays', 'CellLine', 'sampleid'
GENEINFO_CSV = "/homes/hsuleman/pgx_github/data_formatting/gene_info.csv"                  # Contains mappings between 'gene_symbol', 'ensembl_id', 'entrez_id'
OUT_TSV = "/homes/hsuleman/PharmacoGx/analysis/crown_csa_data/pgx_data/pgx_cnv_final.tsv"

# ====== Step 1: Load expression data ======
expr_df = pd.read_csv(EXPR_CSV)
expr_df = expr_df.rename(columns={"Feature": "gene_symbol"})

# ====== Step 2: Map SNP.arrays to CellLine ======
cell_info = pd.read_csv(CELLINFO_CSV)
mapping = dict(zip(cell_info["SNP.arrays"], cell_info["CellLine"]))

# Keep only columns that have a mapping
expr_df = expr_df.loc[:, ["gene_symbol"] + [c for c in expr_df.columns[1:] if c in mapping]]
expr_df = expr_df.rename(columns=mapping)

# ====== Step 3: Transpose so samples are rows, genes are columns ======
expr_df = expr_df.set_index("gene_symbol").T
expr_df.index.name = "Sample"

# ====== Step 4: Convert log2(x) → log2(x+1) ======
expr_df = expr_df.astype(float)
expr_df = np.logaddexp2(0, expr_df)  

# ====== Step 5: Map IDs from gene_info ======
gene_info = pd.read_csv(GENEINFO_CSV)
entrez_map = gene_info.drop_duplicates(subset=["gene_symbol"]).set_index("gene_symbol")["entrez_id"].to_dict()
ensembl_map = gene_info.drop_duplicates(subset=["gene_symbol"]).set_index("gene_symbol")["ensembl_id"].to_dict()

gene_symbols = expr_df.columns
entrez_ids = [entrez_map.get(g, "") for g in gene_symbols]
ensembl_ids = [ensembl_map.get(g, "") for g in gene_symbols]

# ====== Step 6: Drop genes missing any ID ======
keep_mask = [
    str(e).strip() != "" and str(s).strip() != "" and str(en).strip() != ""
    for e, s, en in zip(entrez_ids, gene_symbols, ensembl_ids)
]
expr_df = expr_df.loc[:, keep_mask]
entrez_ids = [e for e, k in zip(entrez_ids, keep_mask) if k]
gene_symbols = [s for s, k in zip(gene_symbols, keep_mask) if k]
ensembl_ids = [en for en, k in zip(ensembl_ids, keep_mask) if k]

# ====== Step 7: Filter out rows not in valid sample IDs ======
valid_sample_ids = cell_info["sampleid"].dropna().astype(str).str.strip().unique()
expr_df.index = expr_df.index.astype(str).str.strip()
expr_df = expr_df.loc[expr_df.index.isin(valid_sample_ids)]

# ====== Step 8: Build final dataframe with gene metadata rows ======
header_rows = pd.DataFrame([
    [""] + entrez_ids,
    [""] + gene_symbols,
    [""] + ensembl_ids
])
data_rows = expr_df.reset_index()
data_rows.columns = ["SampleID"] + gene_symbols

# Match column names for concat
header_rows.columns = data_rows.columns

final_df = pd.concat([header_rows, data_rows], ignore_index=True)


# fill NaNs with 0
final_df = final_df.fillna(0)


# ====== Step 9: Save to TSV ======
final_df.to_csv(OUT_TSV, sep="\t", index=False, header=False)
print(f"Saved final cleaned file to: {OUT_TSV}")
print(f"Genes kept: {len(gene_symbols)}, Samples kept: {expr_df.shape[0]}")
