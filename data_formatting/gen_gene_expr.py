import pandas as pd
import numpy as np
from pathlib import Path

# -----------------------------
# File paths
# -----------------------------
input_folder = Path("/homes/hsuleman/PharmacoSet_Exports3/gene_expr")
gene_info_csv = Path("/homes/hsuleman/PharmacoGx/analysis/outputs/gene_info.csv")

# -----------------------------
# Load gene info mapping
# -----------------------------
gi = pd.read_csv(gene_info_csv, dtype=str)
gi["ensembl_base"] = gi["ensembl_id"].astype(str).str.split(".").str[0]
gi = gi.dropna(subset=["ensembl_base", "gene_symbol", "entrez_id"])
gi = gi[(gi["ensembl_base"] != "") & (gi["gene_symbol"] != "") & (gi["entrez_id"] != "")]
gi = gi.drop_duplicates(subset=["ensembl_base"], keep="first")
mapping = gi.set_index("ensembl_base")[["entrez_id", "gene_symbol"]]

# -----------------------------
# Process each dataset in folder
# -----------------------------
for input_csv in sorted(input_folder.glob("*_rnaseq_summarized.csv")):
    dataset_name = input_csv.name.replace("_rnaseq_summarized.csv", "")
    output_tsv = input_folder / f"{dataset_name}_processed_rnaseq.tsv"

    print(f"Processing {dataset_name}...")

    # 1) Read input
    df = pd.read_csv(input_csv, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")

    # 2) Strip version number ".XX" from Ensemble ID
    df.index = df.index.astype(str).str.split(".").str[0]

    # 4) Convert log2(x+0.001) -> log2(x+1)
    df = pd.DataFrame(
        np.log2(np.power(2.0, df.values) + 0.999),
        index=df.index,
        columns=df.columns
    )

    # 3) Transpose
    expr = df.T

    # 5) Keep only mapped genes
    keep_cols_raw = [g for g in expr.columns if g in mapping.index]
    expr = expr.loc[:, keep_cols_raw]

    # Collapse duplicates
    expr = expr.T.groupby(level=0).mean().T

    # Header rows
    ordered_genes = list(expr.columns)
    cols_out = [""] + ordered_genes
    ensembl_row = pd.DataFrame([[""] + ordered_genes], columns=cols_out)
    entrez_row  = pd.DataFrame([[""] + [mapping.loc[g, "entrez_id"] for g in ordered_genes]], columns=cols_out)
    symbol_row  = pd.DataFrame([[""] + [mapping.loc[g, "gene_symbol"] for g in ordered_genes]], columns=cols_out)

    # Add sample ID as first column
    expr_out = expr.copy()
    expr_out.insert(0, "", expr_out.index)
    expr_out.columns = cols_out

    # Combine and save
    out = pd.concat([ensembl_row, entrez_row, symbol_row, expr_out.reset_index(drop=True)],
                    axis=0, ignore_index=True, sort=False)
    out.to_csv(output_tsv, sep="\t", header=False, index=False)

    print(f"  ✔ Saved to {output_tsv}")

print("All datasets processed.")

