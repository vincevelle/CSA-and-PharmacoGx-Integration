import argparse
import random
from pathlib import Path

import pandas as pd
import numpy as np


TEST_SOURCES = {
    "CCLE_2015",
    "CTRPv2_2015",
    "GDSC_2020(v2-8.2)",
}


def read_list(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def chunk_list(lst, n_chunks):
    """
    Split list into n_chunks as evenly as possible
    """
    return np.array_split(lst, n_chunks)


def write_indices(path, indices):
    with open(path, "w") as f:
        for idx in sorted(indices):
            f.write(f"{idx}\n")


def write_list(path, items):
    with open(path, "w") as f:
        for item in sorted(items):
            f.write(f"{item}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--drugs", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--num_folds", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)

    # -----------------------------
    # Load data
    # -----------------------------
    df = pd.read_csv(args.tsv, sep="\t")
    df.reset_index(drop=True, inplace=True)

    drugs = read_list(args.drugs)
    samples = read_list(args.samples)

    # -----------------------------
    # Create true CV folds
    # -----------------------------
    rng.shuffle(drugs)
    rng.shuffle(samples)

    drug_folds = chunk_list(drugs, args.num_folds)
    sample_folds = chunk_list(samples, args.num_folds)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------
    # Create splits
    # -----------------------------
    for fold_idx in range(args.num_folds):
        split_dir = out_dir / f"split_{fold_idx + 1}"
        split_dir.mkdir(exist_ok=True)

        holdout_drugs = set(drug_folds[fold_idx])
        holdout_samples = set(sample_folds[fold_idx])

        # Drug-blind
        drug_blind = set(
            df.index[df["improve_chem_id"].isin(holdout_drugs)]
        )

        # Sample-blind
        sample_blind = set(
            df.index[df["improve_sample_id"].isin(holdout_samples)]
        )

        blind_indices = drug_blind | sample_blind

        # -----------------------------
        # General test (both-seen)
        # -----------------------------
        eligible_general = df.index[
            (~df.index.isin(blind_indices))
            & (df["source"].isin(TEST_SOURCES))
        ].tolist()

        rng.shuffle(eligible_general)
        n_general = int(0.10 * len(eligible_general))
        general_test = set(eligible_general[:n_general])

        # -----------------------------
        # Train / validation
        # -----------------------------
        remaining = list(
            set(df.index) - blind_indices - general_test
        )
        rng.shuffle(remaining)

        n_val = int(0.10 * len(remaining))
        val_indices = set(remaining[:n_val])
        train_indices = set(remaining[n_val:])

        test_indices = blind_indices | general_test

        # -----------------------------
        # Write outputs
        # -----------------------------
        write_indices(split_dir / "test.txt", test_indices)
        write_indices(split_dir / "train.txt", train_indices)
        write_indices(split_dir / "val.txt", val_indices)
        write_indices(split_dir / "both_seen.txt", general_test)

        write_list(split_dir / "holdout_drugs.txt", holdout_drugs)
        write_list(split_dir / "holdout_samples.txt", holdout_samples)


if __name__ == "__main__":
    main()
