import configparser
import subprocess

# Load param template
param_file = "uno_default_model.txt"  # replace with your param filename
config = configparser.ConfigParser()
config.read(param_file)

for fold in range(1, 11):  # 1 through 10
    # Update split files
    config["Preprocess"]["train_split_file"] = f"example_splits/split_{fold}/train.txt"
    config["Preprocess"]["val_split_file"] = f"example_splits/split_{fold}/val.txt"
    config["Preprocess"]["test_split_file"] = f"example_splits/split_{fold}/test.txt"
    config["Preprocess"]["ml_data_outdir"] = f"results/fold_{fold}"

    # Save back to params.ini (overwrites)
    with open(param_file, "w") as f:
        config.write(f)

    # Build the command with fold-specific dirs
    outdir = f"results/fold_{fold}"
    cmd = (
        f"python uno_preprocess_improve.py --input_dir new_benchmark_dataset --output_dir {outdir} "
        f"&& python uno_train_improve.py --input_dir {outdir} --output_dir {outdir} "
        f"&& python uno_infer_improve.py --input_data_dir {outdir} --input_model_dir {outdir} "
        f"--output_dir {outdir} --calc_infer_score true"
    )
    
    print(f"\n=== Running Fold {fold} ===")
    subprocess.run(cmd, shell=True, check=True)
