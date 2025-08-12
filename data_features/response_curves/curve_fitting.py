import os
import pandas as pd
import numpy as np
import time
from multiprocessing import Pool
from tqdm import tqdm
from drug_response_curve import drug_screening_dataframe

# === Settings ===
DATA_DIR = "/homes/hsuleman/PharmacoGx/PharmacoSet_Exports"   # Folder with PharmacoGx data saved as CSVs
OUT_DIR = "./dose_response_output_FULL"

SELECTED_DATASETS = [
    "NCI60_2021",
    "PRISM_2020",
    "GDSC_2020(v1-8.2)",
    "GDSC_2020(v2-8.2)",
    "CTRPv2_2015",
    "CCLE_2015",
    "FIMM_2016",
    "gCSI_2019"
]

NUM_WORKERS = 65   # Initial run was on an 80 CPU machine, adjust according to your machine spec.
SKIP_FIRST_N = 0 
PLOT = False

# param_bound = ([low_Einf, low_EC50, low_HS, low_E0], [up_Einf, up_EC50, up_HS, up_E0])
PARAM_BOUND = ([0, -19, -5, 1], [1, 5, 5, 1 + 1e-6])  

TIMING_LOG = os.path.join(OUT_DIR, "dataset_timings.tsv") # Keeps track of how long each dataset takes to generate

header = [
    "source", "sample_id", "drug_id", "exp_id", "auc", "ic50", "ec50", "ec50se",
    "r2fit", "einf", "hs", "aac1", "auc1", "flipped", "R2_orig", "R2_flipped"
]

def load_long_format(dataset_name, skip_first_n=SKIP_FIRST_N):
    file_path = os.path.join(DATA_DIR, dataset_name, f"{dataset_name}_sensitivity_raw.csv")
    print(f"Loading dataset: {dataset_name}")
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return pd.DataFrame()

    try:
        t0 = time.time()
        df = pd.read_csv(file_path, index_col=0, low_memory=False)
        print(f"CSV loaded in {time.time() - t0:.2f}s | Shape: {df.shape}")
    except Exception as e:
        print(f"Failed to load {dataset_name}: {e}")
        return pd.DataFrame()

    # Skip the first N unique experiments
    all_exp_ids = df.index.unique()
    if len(all_exp_ids) <= skip_first_n:
        print(f"Dataset only has {len(all_exp_ids)} experiments; skipping.")
        return pd.DataFrame()

    selected_ids = all_exp_ids[skip_first_n:]
    df = df.loc[selected_ids]
    print(f"Selected {len(selected_ids)} experiments after skipping first {skip_first_n}")

    dose_cols = [col for col in df.columns if "Dose" in col and "dose." in col]
    all_rows = []
    t1 = time.time()
    for i in range(1, len(dose_cols) + 1):
        try:
            temp = df[['treatmentid', 'sampleid']].copy()
            temp['dose'] = df[f'dose.{i}.Dose']
            temp['viability'] = df[f'dose.{i}.Viability']
            temp['dose_log10M'] = np.log10(temp['dose'].astype(float) * 1e-6)  # PharmacoGx doses are initially in micromolar, so we convert to molar
            temp['viability_norm'] = np.clip(temp['viability'], 0, 200) / 100.0  # Truncates viabilities to range (0,200), since some datasets had extreme outliers
            temp['exp.id'] = df.index   # Adjust to name of column, if first col of sensitivity file has named 'Experiment ID'
            temp['plate_id'] = 'P1'
            temp['is_control'] = False
            temp['dataset'] = dataset_name
            all_rows.append(temp)
        except Exception as e:
            print(f"Skipping dose.{i} due to error: {e}")
    long_df = pd.concat(all_rows, ignore_index=True)
    print(f"Long-format conversion done in {time.time() - t1:.2f}s | Final shape: {long_df.shape}")
    return long_df

def process_experiment(args):
    exp_id, dataset_name, df_exp, plot_dir = args
    try:
        dsd = drug_screening_dataframe(
            res=df_exp,
            control_types=[],
            model_col='sampleid',
            plate_col='plate_id',
            row_col=None,
            column_col=None,
            comp_col='treatmentid',
            dose_col='dose_log10M',
            solvent_col=None,
            fea_col='viability',
            across_plate_col=None,
            study_id=dataset_name,
            exp_col='exp.id',
            norm_fea_col='viability_norm',
            control_col='is_control'
        )

        results, _ = dsd.fit_dose_response_curve(param_bound=PARAM_BOUND, plot_save_folder=plot_dir if PLOT else "./")

        rows = []
        for key, fit_obj in results.items():
            metrics = fit_obj.fit_metrics.copy()
            rows.append({
                "source": fit_obj.study_id,
                "sample_id": fit_obj.specimen,
                "drug_id": fit_obj.compound,
                "exp_id": fit_obj.exp_id,
                "auc": metrics.get("AUC"),
                "ic50": metrics.get("IC50"),
                "ec50": metrics.get("EC50"),
                "ec50se": metrics.get("EC50_se"),
                "r2fit": metrics.get("R2fit"),
                "einf": metrics.get("Einf"),
                "hs": metrics.get("HS"),
                "aac1": metrics.get("AAC1"),
                "auc1": metrics.get("AUC1"),
                "flipped": metrics.get("flipped", False),   #  boolean flag to check if curve was flipped
                "R2_orig": metrics.get("R2_orig", np.nan),    #  R^2 without flipping
                "R2_flipped": metrics.get("R2_flipped", np.nan)   # R^2 with curve flipping
            })
        return rows
    except Exception as e:
        print(f"Curve fitting failed for experiment {exp_id} in {dataset_name}: {e}")
        return []

def process_dataset(ds):
    print(f"\nStarting dataset: {ds}")
    start_time = time.time()

    df_long = load_long_format(ds)
    if df_long.empty:
        print(f"Skipping empty dataset: {ds}")
        return

    exp_ids = df_long['exp.id'].unique()
    print(f"Unique experiments in long data: {len(exp_ids)}")

    # Efficient pre-slicing using groupby, prevents copying large datasets between processes
    print(f"Grouping by exp.id ...")
    t_group = time.time()
    grouped_experiments = dict(tuple(df_long.groupby('exp.id')))
    print(f"Grouped in {time.time() - t_group:.2f}s")

    plot_dir = os.path.join(OUT_DIR, ds, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    output_file = os.path.join(OUT_DIR, f"{ds}_responseFULL.tsv")
    if os.path.exists(output_file):
        print(f"Skipping {ds} (already processed)")
        return

    task_list = [
        (exp_id, ds, grouped_experiments[exp_id], plot_dir)
        for exp_id in exp_ids if exp_id in grouped_experiments
    ]

    print(f"Starting multiprocessing for {len(task_list)} experiments...")
    with open(output_file, 'w') as f_out:
        f_out.write('\t'.join(header) + '\n')
        with Pool(NUM_WORKERS) as pool:
            for res in tqdm(pool.imap_unordered(process_experiment, task_list, chunksize=10),
                            total=len(task_list), desc=f"{ds} (full run)"):
                if res:
                    pd.DataFrame(res).to_csv(f_out, sep='\t', index=False, header=False)

    duration = time.time() - start_time
    print(f"{ds} full run completed in {duration:.2f} seconds")

    with open(TIMING_LOG, 'a') as log:
        log.write(f"{ds}\t{len(task_list)}\t{duration:.2f}\n")

def init_timing_log():
    os.makedirs(OUT_DIR, exist_ok=True)
    if not os.path.exists(TIMING_LOG):
        with open(TIMING_LOG, 'w') as f:
            f.write("dataset\texperiments\tduration_sec\n")

if __name__ == "__main__":
    init_timing_log()
    for ds in SELECTED_DATASETS:
        process_dataset(ds)
