import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error

# === CONFIG ===
file1 = "/homes/hsuleman/PharmacoGx/analysis/scripts/labeled_response.tsv"
file2 = "/homes/hsuleman/PharmacoGx/analysis/scripts/dose_response_output_FULL/final_responses/pgx_response.tsv"
output_matches = "overlapping_experiments.tsv"

source_map = {
    "CCLE": "CCLE_2015",
    "CTRPv2": "CTPRv2_2015",
    "gCSI": "gCSI_2019",
    "GDSCv1": "GDSC_2020(v1-8.2)",
    "GDSCv2": "GDSC_2020(v2-8.2)"
}

def normalize_source(source):
    return source_map.get(source, source)

# === LOAD AND PREPROCESS ===
df1 = pd.read_csv(file1, sep="\t", dtype=str)
df2 = pd.read_csv(file2, sep="\t", dtype=str)

df2 = df2.rename(columns={"sample_id": "improve_sample_id", "drug_id": "improve_chem_id"})

df1["source_norm"] = df1["source"].map(normalize_source)
df2["source_norm"] = df2["source"].map(normalize_source)

df1["match_key"] = df1["source_norm"] + "|" + df1["improve_sample_id"] + "|" + df1["improve_chem_id"]
df2["match_key"] = df2["source_norm"] + "|" + df2["improve_sample_id"] + "|" + df2["improve_chem_id"]

# === FIND OVERLAPS ===
overlapping_keys = set(df1["match_key"]) & set(df2["match_key"])

# Deduplicate per file
df1_dedup = (
    df1[df1["match_key"].isin(overlapping_keys)]
    .drop_duplicates(subset="match_key")
    .copy()
)
df2_dedup = (
    df2[df2["match_key"].isin(overlapping_keys)]
    .drop_duplicates(subset="match_key")
    .copy()
)

# Convert AUC columns to numeric
df1_dedup["auc_file1"] = pd.to_numeric(df1_dedup["auc"], errors="coerce")
df2_dedup["auc_file2"] = pd.to_numeric(df2_dedup["auc"], errors="coerce")

# === Merge on match_key ===
merged = pd.merge(
    df1_dedup[["match_key", "source_norm", "improve_sample_id", "improve_chem_id", "auc_file1"]],
    df2_dedup[["match_key", "auc_file2"]],
    on="match_key",
    how="inner"
)

# Drop invalid AUCs
merged = merged.replace([np.inf, -np.inf], np.nan).dropna(subset=["auc_file1", "auc_file2"])

# === SAVE OUTPUT FILE ===
output_df = merged[["source_norm", "improve_sample_id", "improve_chem_id", "auc_file1", "auc_file2"]]
output_df = output_df.rename(columns={"source_norm": "source"})
output_df.to_csv(output_matches, sep="\t", index=False)

# === CALCULATE AUC METRICS ===
x_valid = output_df["auc_file1"]
y_valid = output_df["auc_file2"]

print("Correlation and error metrics for AUC (based on overlapping experiments):\n")

if len(output_df) < 2:
    print("Not enough valid AUC data points to compute metrics.")
else:
    pearson, _ = pearsonr(x_valid, y_valid)
    spearman, _ = spearmanr(x_valid, y_valid)
    mae = mean_absolute_error(x_valid, y_valid)
    rmse = np.sqrt(mean_squared_error(x_valid, y_valid))

    print(f"AUC:")
    print(f"  Pearson  = {pearson:.4f}")
    print(f"  Spearman = {spearman:.4f}")
    print(f"  MAE      = {mae:.4f}")
    print(f"  RMSE     = {rmse:.4f}")
