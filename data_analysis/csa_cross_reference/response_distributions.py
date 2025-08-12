import pandas as pd
import matplotlib.pyplot as plt

# Load TSV file
file_path = '/homes/hsuleman/PharmacoGx/analysis/scripts/dose_response_output_FULL/final_responses/pgx_response.tsv'  # <-- Replace with file to analyze
df = pd.read_csv(file_path, sep='\t')

# Summary statistics for all AUC values
auc_all = df['auc'].dropna()
print("=== Summary Statistics for ALL AUC values ===")
print(auc_all.describe())
print(f"Standard Deviation: {auc_all.std():.4f}")
print()

# Filter AUCs between 0 and 1
auc_0_1 = auc_all[(auc_all > 0) & (auc_all < 1)]
print("=== Summary Statistics for AUC values between 0 and 1 ===")
print(auc_0_1.describe())
print(f"Standard Deviation: {auc_0_1.std():.4f}")
print()

# Count and percentage of out-of-bounds values
total_count = len(auc_all)
greater_than_1 = (auc_all > 1).sum()
less_than_0 = (auc_all < 0).sum()
greater_than_10 = (auc_all > 10).sum()

print("=== Out-of-Range AUC Values ===")
print(f">1: {greater_than_1} ({greater_than_1 / total_count * 100:.2f}%)")
print(f"<0: {less_than_0} ({less_than_0 / total_count * 100:.2f}%)")
print(f">10: {greater_than_10} ({greater_than_10 / total_count * 100:.2f}%)")
print()

# Histogram for AUC values in the range (0, 1)
plt.figure(figsize=(10, 6))
plt.hist(auc_all, bins=50, edgecolor='black')  #
plt.title('PharmacoGx Response')
plt.xlabel('AUC')
plt.ylabel('Frequency')
plt.grid(True)
plt.tight_layout()

histogram_path = '/homes/hsuleman/PharmacoGx/analysis/auc__full_corrected2'
plt.savefig(histogram_path)
print(f"Histogram saved as: {histogram_path}")

plt.show()