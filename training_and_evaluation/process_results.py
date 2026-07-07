import os
import json
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr, spearmanr


# -----------------------------
# CONFIG
# -----------------------------
SPLITS_DIR = Path("results") # add path to results
REFERENCE_TSV = Path("response.tsv") # add path to response file

ORGANOID_SOURCES = {
    "Lee",
    "Narasimhan",
    "Tempus",
    "VanDerWeltering"
}


# -----------------------------
# Utility functions
# -----------------------------
def read_list(path):
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())


def compute_metrics(df):
    y_true = df["auc_true"].values
    y_pred = df["auc_pred"].values

    mse = mean_squared_error(y_true, y_pred)
    return {
        "num_samples": int(len(df)),
        "MSE": float(mse),
        "RMSE": float(np.sqrt(mse)),
        "PCC": float(pearsonr(y_true, y_pred)[0]),
        "SCC": float(spearmanr(y_true, y_pred)[0]),
        "R2": float(r2_score(y_true, y_pred)),
    }


# -----------------------------
# Load reference once
# -----------------------------
ref = pd.read_csv(REFERENCE_TSV, sep="\t")
ref["auc"] = ref["auc"].round(4)

ref_grouped = (
    ref[["source", "improve_sample_id", "improve_chem_id", "auc"]]
    .groupby(["improve_sample_id", "improve_chem_id", "auc"], as_index=False)
    .agg({"source": lambda x: sorted(set(x))})
)


# -----------------------------
# Main loop over folds
# -----------------------------
all_fold_metrics = {}

for fold_dir in sorted(SPLITS_DIR.glob("fold_*")):

    fold_num = int(fold_dir.name.split("_")[1])

    print(f"\n=== Processing {fold_dir.name} ===")

    final_eval_dir = fold_dir / "final_eval"
    final_eval_dir.mkdir(exist_ok=True)

    # -----------------------------
    # Load predictions
    # -----------------------------
    pred_path = fold_dir / "test_y_data_predicted_filtered.csv"
    df = pd.read_csv(pred_path)
    df["auc_true"] = df["auc_true"].round(4)

    # -----------------------------
    # Attach source info
    # -----------------------------
    # df = df.merge(
    #     ref_grouped,
    #     on=["improve_sample_id", "improve_chem_id", "auc"],
    #     how="left"
    # )
    # df["source"] = df["source"].apply(
    #     lambda x: ";".join(x) if isinstance(x, list) else np.nan
    # )

    # # -----------------------------
    # # Expand compounded source
    # # -----------------------------
    # df["source"] = df["source"].str.split(";")
    # df = df.explode("source")
    # df["source"] = df["source"].str.strip()

    # -----------------------------
    # Organoid split
    # -----------------------------
    df_org = df[df["source"].isin(ORGANOID_SOURCES)]
    df_non_org = df[~df["source"].isin(ORGANOID_SOURCES)]

    df_org.to_csv(final_eval_dir / "ORGANOIDS.csv", index=False)
    df_non_org.to_csv(final_eval_dir / "NO_ORGANOIDS.csv", index=False)

    # -----------------------------
    # Load holdouts
    # -----------------------------
    holdout_drugs = read_list(fold_dir / "holdout_drugs.txt")
    holdout_samples = read_list(fold_dir / "holdout_samples.txt")

    df_non_org["is_holdout_drug"] = df_non_org["improve_chem_id"].isin(holdout_drugs)
    df_non_org["is_holdout_sample"] = df_non_org["improve_sample_id"].isin(holdout_samples)

    # -----------------------------
    # Split by evaluation strategy
    # -----------------------------
    splits = {
        "general": df_non_org[(~df_non_org.is_holdout_drug) & (~df_non_org.is_holdout_sample)],
        "holdout_sample": df_non_org[(~df_non_org.is_holdout_drug) & (df_non_org.is_holdout_sample)],
        "holdout_drug": df_non_org[(df_non_org.is_holdout_drug) & (~df_non_org.is_holdout_sample)],
        "disjoint": df_non_org[(df_non_org.is_holdout_drug) & (df_non_org.is_holdout_sample)],
    }

    fold_metrics = {}

    for name, subdf in splits.items():
        out_csv = final_eval_dir / f"{name}_predictions.csv"
        subdf.drop(columns=["is_holdout_drug", "is_holdout_sample"]).to_csv(out_csv, index=False)

        if len(subdf) > 1:
            fold_metrics[name] = compute_metrics(subdf)

    # -----------------------------
    # Save per-fold metrics
    # -----------------------------
    metrics_path = final_eval_dir / "metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(fold_metrics, f, indent=2)

    all_fold_metrics[fold_dir.name] = fold_metrics


# -----------------------------
# Save global summary
# -----------------------------
summary_path = SPLITS_DIR / "all_folds_metrics.json"
with open(summary_path, "w") as f:
    json.dump(all_fold_metrics, f, indent=2)

print("\n All folds processed successfully.")
