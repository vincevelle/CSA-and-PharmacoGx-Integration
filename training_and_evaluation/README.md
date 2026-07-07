# Cross-Validation Pipeline

This folder contains helper scripts for creating 10-fold cross-validation splits, training models, processing fold-level predictions, and plotting aggregate performance across evaluation settings.

See https://github.com/JDACS4C-IMPROVE/UNO and https://github.com/JDACS4C-IMPROVE/GraphDRP/tree/develop for instructions on how to run the models used in this analysis

The pipeline is organized around four scripts:

| Script | Purpose |
|---|---|
| `gen_splits.py` | Creates cross-validation split files for drug-blind, sample-blind, disjoint, and mixed set testing. |
| `run_cv.py` | Runs preprocessing, training, and inference for each fold using the generated split files. |
| `process_results.py` | Separates predictions into evaluation categories and computes metrics for each fold. |
| `radar_plots.py` | Generates a radar plot comparing model performance across evaluation categories. |

---

## Overview

The benchmark evaluates model generalization under four settings:

1. **Mixed**  
   The drug and cancer sample have both appeared in training, but the specific pair is held out.

2. **Drug-Blind**  
   The drug is held out from training, while the cancer sample may be seen.

3. **Cancer-Blind**  
   The cancer sample is held out from training, while the drug may be seen.

4. **Disjoint**  
   Both the drug and cancer sample are held out from training.

The split-generation script creates fold-specific holdout drugs and holdout samples. Rows involving either held-out drugs or held-out samples are placed into the test set. 10% of remaining rows from selected test sources are also added to the test set.

---

## Expected Input Files

### Response TSV

`gen_splits.py` expects a tab-separated response file containing at least the following columns:

```text
improve_sample_id
improve_chem_id
source
auc
```

Example:

```text
improve_sample_id	improve_chem_id	source	auc
ACH-000001	Drug_001	CCLE_2015	0.7345
ACH-000002	Drug_002	GDSC_2020(v2-8.2)	0.5123
```

### Drug list

A plain text file containing one drug ID per line:

```text
PX-12
Triptolide
Prima-1
AM-580
```

### Sample list

A plain text file containing one sample ID per line:

```text
CVCL_1086
CVCL_2010
CVCL_2791
CVCL_2420
```

---



## Step 1: Generate Cross-Validation Splits

Run:

```bash
python gen_splits.py \
  --tsv response.tsv \
  --drugs drugs.txt \
  --samples samples.txt \
  --out_dir example_splits \
  --num_folds 10 \
  --seed 42
```

This creates:

```text
example_splits/split_1/
example_splits/split_2/
...
example_splits/split_10/
```

Each split directory contains:

| File | Description |
|---|---|
| `train.txt` | Row indices used for training. |
| `val.txt` | Row indices used for validation. |
| `test.txt` | Row indices used for testing. |
| `both_seen.txt` | Mixed set test rows where both drug and sample are seen in training, but not their pairing. |
| `holdout_drugs.txt` | Drug IDs excluded from training for that fold. |
| `holdout_samples.txt` | Sample IDs excluded from training for that fold. |

### Test sources used for the general test set

The script samples test rows from:

```python
TEST_SOURCES = {
    "CCLE_2015",
    "CTRPv2_2015",
    "GDSC_2020(v2-8.2)",
}
```

Rows from these sources are eligible for the mixed test subset if they are not already part of a drug-blind or sample-blind holdout.

---

## Step 2: Run Cross-Validation

The script expects generated splits at:

```text
example_splits/split_1/
example_splits/split_2/
...
example_splits/split_10/
```

It writes fold outputs to:

```text
results/fold_1/
results/fold_2/
...
results/fold_10/
```

Run:

```bash
python run_cv.py
```

For each fold, the script runs:

```bash
python uno_preprocess_improve.py --input_dir new_benchmark_dataset --output_dir results/fold_N
python uno_train_improve.py --input_dir results/fold_N --output_dir results/fold_N
python uno_infer_improve.py --input_data_dir results/fold_N --input_model_dir results/fold_N --output_dir results/fold_N --calc_infer_score true
```

where `N` is the fold number.

---

## Step 3: Copy Holdout Metadata Into Fold Results

`process_results.py` expects these files to exist inside each `results/fold_N/` directory:

```text
holdout_drugs.txt
holdout_samples.txt
```

However, `gen_splits.py` writes them to `example_splits/split_N/`.

Before processing results, copy them into the corresponding result folders:

```bash
for i in {1..10}; do
  cp example_splits/split_${i}/holdout_drugs.txt results/fold_${i}/holdout_drugs.txt
  cp example_splits/split_${i}/holdout_samples.txt results/fold_${i}/holdout_samples.txt
done
```

This step is required so that `process_results.py` can identify whether each prediction belongs to the mixed, drug-blind, sample-blind, or disjoint evaluation category.

---

## Step 4: Process Fold Results

Edit the configuration section in `process_results.py` if needed:

```python
SPLITS_DIR = Path("results")
REFERENCE_TSV = Path("response.tsv")
```

The script expects each fold directory to contain a prediction file named:

```text
test_y_data_predicted.csv
```

The prediction CSV should contain at least:

```text
improve_sample_id
improve_chem_id
auc_true
auc_pred
source
```

Run:

```bash
python process_results.py
```

For each fold, this creates:

```text
results/fold_N/final_eval/
├── ORGANOIDS.csv
├── NO_ORGANOIDS.csv
├── general_predictions.csv
├── holdout_sample_predictions.csv
├── holdout_drug_predictions.csv
├── disjoint_predictions.csv
└── metrics.json
```

It also writes a global metrics summary:

```text
results/all_folds_metrics.json
```

---

## Evaluation Categories

After organoid rows are removed, `process_results.py` assigns non-organoid predictions to one of four groups:

| Output file | Meaning |
|---|---|
| `general_predictions.csv` | Neither drug nor sample is held out. |
| `holdout_sample_predictions.csv` | Sample is held out, drug is not held out. |
| `holdout_drug_predictions.csv` | Drug is held out, sample is not held out. |
| `disjoint_predictions.csv` | Both drug and sample are held out. |

Organoid rows are identified using:

```python
ORGANOID_SOURCES = {
    "Lee",
    "Narasimhan",
    "Tempus",
    "VanDerWeltering"
}
```

These rows are written separately to `ORGANOIDS.csv` and excluded from the four standard non-organoid evaluation categories.

---

## Metrics

For each evaluation category, the following metrics are computed:

| Metric | Description |
|---|---|
| `num_samples` | Number of evaluated rows. |
| `MSE` | Mean squared error. |
| `RMSE` | Root mean squared error. |
| `PCC` | Pearson correlation coefficient. |
| `SCC` | Spearman correlation coefficient. |
| `R2` | Coefficient of determination. |



---

## Step 5: Generate Radar Plot

After computing fold-level results, update the values in `radar_plots.py`:

```python
gdrp_old = [.75, -.11, .60, -.26]
gdrp_new = [.66, .12, .57, .08]
```

Values should be ordered as:

```python
["Mixed", "Drug-Blind", "Cancer-Blind", "Disjoint"]
```

Run:

```bash
python radar_plots.py
```

The script writes:

```text
gdrp_radarplot.png
```

---

## Important Notes

### Source column in predictions

`process_results.py` currently assumes the prediction CSV already contains a `source` column.

If your prediction file does not include `source`, uncomment the merge block in `process_results.py` so source information is recovered from `response.tsv`.

### AUC rounding

Both reference AUC and prediction true AUC are rounded to four decimals in `process_results.py`:

```python
ref["auc"] = ref["auc"].round(4)
df["auc_true"] = df["auc_true"].round(4)
```

This helps align prediction rows with reference rows when source information needs to be merged.

### Reproducibility

Splits are randomized using the seed provided to `gen_splits.py`:

```bash
--seed 42
```

Use the same seed to regenerate identical split assignments.

---

## Full Pipeline Example

```bash
# 1. Generate split files
python gen_splits.py \
  --tsv response.tsv \
  --drugs drugs.txt \
  --samples samples.txt \
  --out_dir example_splits \
  --num_folds 10 \
  --seed 42

# 2. Run preprocessing, training, and inference for each fold
python run_cv.py

# 3. Copy holdout metadata into fold result directories
for i in {1..10}; do
  cp example_splits/split_${i}/holdout_drugs.txt results/fold_${i}/holdout_drugs.txt
  cp example_splits/split_${i}/holdout_samples.txt results/fold_${i}/holdout_samples.txt
done

# 4. Process predictions and compute metrics
python process_results.py

# 5. Generate radar plot after updating metric values
python radar_plots.py
```

---

## Outputs

After the full pipeline completes, the main outputs are:

```text
results/all_folds_metrics.json
results/fold_*/final_eval/metrics.json
results/fold_*/final_eval/*_predictions.csv
gdrp_radarplot.png
```

These files can be used to compare model performance across mixed, drug-blind, cancer-blind, and disjoint generalization settings.
