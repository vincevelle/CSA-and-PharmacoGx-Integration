import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json

# Base directory
base_dir = "PharmacoSet_Exports"

# Output directory
output_dir = os.path.join(base_dir, "histogram_outputs")
os.makedirs(output_dir, exist_ok=True)

# Fixed histogram ranges
log_dose_range = (-3.5, 2)
viability_range = (0, 200)

# Binning
bins_dose = np.linspace(*log_dose_range, 50)
bins_viability = np.linspace(*viability_range, 50)

# Outlier logs
dose_outliers = []
viability_outliers = []

# Loop through each dataset folder
for folder in os.listdir(base_dir):
    folder_path = os.path.join(base_dir, folder)
    if not os.path.isdir(folder_path):
        continue

    dataset_name = folder
    csv_path = os.path.join(folder_path, f"{dataset_name}_sensitivity_raw.csv")

    if not os.path.exists(csv_path):
        continue

    print(f"Processing: {dataset_name}")

    try:
        df = pd.read_csv(csv_path)

        # Identify dose and viability columns
        dose_cols = [col for col in df.columns if col.lower().endswith(".dose")]
        viability_cols = [col for col in df.columns if col.lower().endswith(".viability")]

        # Process doses
        dose_values = df[dose_cols].apply(pd.to_numeric, errors='coerce').values.flatten()
        dose_values = dose_values[~np.isnan(dose_values)]
        dose_values = dose_values[dose_values > 0]
        log_doses = np.log10(dose_values)

        out_doses = log_doses[(log_doses < log_dose_range[0]) | (log_doses > log_dose_range[1])]
        if len(out_doses) > 0:
            print(f"[{dataset_name}] Found {len(out_doses)} log10(dose) outliers")
            dose_outliers.append({
                "dataset": dataset_name,
                "outliers": out_doses.tolist()
            })

        # Plot dose histogram
        plt.figure(figsize=(10, 4))
        plt.hist(log_doses, bins=bins_dose)
        plt.title(f"{dataset_name} - log10(Dose) Distribution")
        plt.xlabel("log10(Dose)")
        plt.ylabel("Frequency")
        plt.xlim(log_dose_range)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{dataset_name}_log10_dose_hist.png"))
        plt.close()

        # Process viabilities
        viab_values = df[viability_cols].apply(pd.to_numeric, errors='coerce').values.flatten()
        viab_values = viab_values[~np.isnan(viab_values)]

        out_viab = viab_values[(viab_values < viability_range[0]) | (viab_values > viability_range[1])]
        if len(out_viab) > 0:
            print(f"[{dataset_name}] Found {len(out_viab)} viability outliers")
            viability_outliers.append({
                "dataset": dataset_name,
                "outliers": out_viab.tolist()
            })

        # Plot viability histogram
        plt.figure(figsize=(10, 4))
        plt.hist(viab_values, bins=bins_viability)
        plt.title(f"{dataset_name} - Viability Distribution")
        plt.xlabel("Viability")
        plt.ylabel("Frequency")
        plt.xlim(viability_range)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{dataset_name}_viability_hist.png"))
        plt.close()

    except Exception as e:
        print(f"Error processing {dataset_name}: {e}")

# Save outlier logs
with open(os.path.join(output_dir, "dose_outliers.json"), "w") as f:
    json.dump(dose_outliers, f, indent=2)

with open(os.path.join(output_dir, "viability_outliers.json"), "w") as f:
    json.dump(viability_outliers, f, indent=2)

print("All histograms generated and outliers logged.")
