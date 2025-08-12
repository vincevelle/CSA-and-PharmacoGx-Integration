def discretize(value):
    try:
        value = float(value)
    except ValueError:
        return ''  # Return blank if non-numeric or missing

    if value < 0.5210507:
        return -2  # Deep deletion
    elif value < 0.7311832:
        return -1  # Het loss
    elif value < 1.214125:
        return 0   # Diploid
    elif value < 1.422233:
        return 1   # Gain
    else:
        return 2   # Amplification

def main(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # First three lines are headers
    position_headers = lines[0].strip().split('\t')
    gene_symbols = lines[1].strip()
    gene_ids = lines[2].strip()

    expected_value_cols = len(position_headers) - 1  # exclude sample name

    data_lines = lines[3:]
    discretized_lines = []

    
    for line_num, line in enumerate(data_lines, start=4):
        parts = line.strip().split('\t')
        sample_name = parts[0]
        values = parts[1:]

        if len(values) < expected_value_cols:
            print(f"⚠️ Warning at line {line_num}: {len(values)} values (expected {expected_value_cols})")
            # Fill missing values with '0'
            values += ['0'] * (expected_value_cols - len(values))

        # Discretize values, treating missing values as '0'
        discretized_values = [str(discretize(val)) for val in values]
        discretized_line = '\t'.join([sample_name] + discretized_values)
        discretized_lines.append(discretized_line)



    # Write output
    with open(output_file, 'w') as f:
        f.write('\t'.join(position_headers) + '\n')
        f.write(gene_symbols + '\n')
        f.write(gene_ids + '\n')
        for line in discretized_lines:
            f.write(line + '\n')

if __name__ == "__main__":
    input_path = "/homes/hsuleman/PharmacoGx/analysis/crown_csa_data/pgx_data/pgx_cnv_final.tsv"             
    output_path = "/homes/hsuleman/PharmacoGx/analysis/crown_csa_data/pgx_data/pgx_discretized_cnv_final.tsv"
    main(input_path, output_path)
