import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# === INPUTS ===
RAW_DATA_DIR = "/homes/hsuleman/PharmacoGx/PharmacoSet_Exports"
AUC_FILE = "/homes/hsuleman/PharmacoGx/analysis/crown_csa_data/pgx_data_uno/y_data/response.tsv"  # path to AUC file
OUTPUT_FILE = "avg_viability_for_auc_gt1.tsv"

# === STEP 1: Load AUC results ===
auc_df = pd.read_csv(AUC_FILE, sep="\t")
print(f"Loaded AUC file: {auc_df.shape}")

# Filter for AUC > 1
auc_gt1 = auc_df[auc_df["auc"] > 1].copy()
print(f"Experiments with AUC > 1: {auc_gt1.shape[0]}")

# Extract dataset names
datasets_with_hits = auc_gt1["source"].unique()

# === STEP 2: Function to load and normalize viability from raw file ===
def load_long_format(dataset_name):
    file_path = os.path.join(RAW_DATA_DIR, dataset_name, f"{dataset_name}_fixed_sensitivity_raw.csv")
    if not os.path.exists(file_path):
        print(f"Missing file: {file_path}")
        return pd.DataFrame()

    df = pd.read_csv(file_path, index_col=0, low_memory=False)
    dose_cols = [col for col in df.columns if "Dose" in col and "dose." in col]
    all_rows = []
    for i in range(1, len(dose_cols) + 1):
        try:
            temp = df[['treatmentid', 'sampleid']].copy()
            temp['dose'] = df[f'dose.{i}.Dose']
            temp['viability'] = df[f'dose.{i}.Viability']
            temp['dose_log10M'] = (temp['dose'].astype(float) * 1e-6).apply(
                lambda x: pd.NA if pd.isna(x) else np.log10(x)
            )
            temp['viability_norm'] = temp['viability'].clip(0, 200) / 100.0
            temp['exp.id'] = df.index
            temp['dataset'] = dataset_name
            all_rows.append(temp)
        except Exception as e:
            print(f"Skipping dose.{i} due to: {e}")
    return pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()

# === STEP 3: Iterate over datasets and collect results ===
result_rows = []

for ds in datasets_with_hits:
    print(f"\nProcessing dataset: {ds}")
    df_long = load_long_format(ds)
    if df_long.empty:
        continue

    exp_ids = auc_gt1.loc[auc_gt1["source"] == ds, "study"].unique()
    subset = df_long[df_long["exp.id"].isin(exp_ids)]
    if subset.empty:
        print(f"No matching experiments in raw data for {ds}")
        continue

    # Average viability per experiment
    avg_viab = subset.groupby("exp.id")["viability_norm"].mean().reset_index()
    avg_viab["dataset"] = ds

    # Merge AUC values into this dataset's averages
    merged = avg_viab.merge(
        auc_gt1[auc_gt1["source"] == ds][["study", "auc"]],
        left_on="exp.id",
        right_on="study",
        how="left"
    )
    merged.drop(columns="study", inplace=True)
    result_rows.append(merged)

# === STEP 4: Save combined results ===
if result_rows:
    final_df = pd.concat(result_rows, ignore_index=True)
    final_df.rename(columns={"viability_norm": "avg_viability_norm"}, inplace=True)
    final_df = final_df[["dataset", "exp.id", "avg_viability_norm", "auc"]]
    final_df.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"\nSaved results to {OUTPUT_FILE} | Rows: {final_df.shape[0]}")

    # === STEP 5: Histogram with thresholds ===
    thresholds = [0.1, 0.25, 0.5, 0.75, 1.0]

    plt.figure(figsize=(10, 6))
    counts, bins, patches = plt.hist(
        final_df["avg_viability_norm"].dropna(),
        bins=50,
        range=(0, 2),
        edgecolor='black'
    )

    plt.xlabel("Average Viability (normalized)")
    plt.ylabel("Count")
    plt.title("Histogram of Average Viability for AUC > 1 Experiments")

    # Add threshold lines & labels
    for th in thresholds:
        plt.axvline(x=th, color='red', linestyle='--', linewidth=1)
        plt.text(th, max(counts)*0.95, f"{th}", rotation=90,
                 verticalalignment='top', horizontalalignment='right',
                 fontsize=8, color='red')

    plt.tight_layout()
    plt.savefig("avg_viability_histogram.png", dpi=300)
    plt.close()
    print("Histogram with thresholds saved as avg_viability_histogram.png")


    # === STEP 6: Counts & percentages for thresholds ===
    thresholds = [0.1, 0.25, 0.5, 0.75, 1.0]
    total = len(final_df)
    print("\nCounts and percentages below thresholds:")
    for th in thresholds:
        count = (final_df["avg_viability_norm"] < th).sum()
        perc = (count / total) * 100
        print(f"  < {th:.2f}: {count} ({perc:.2f}%)")
else:
    print("No results to save.")

