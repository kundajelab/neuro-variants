import pandas as pd
import os
import gzip

# Define the absolute paths
input_folder = '/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/trevino_2021/peaks/overlap'
output_folder = '/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/trevino_2021/peaks/dedup_overlap'

# Ensure the output directory exists
os.makedirs(output_folder, exist_ok=True)

# Process each .bed.gz file in the input directory
for filename in os.listdir(input_folder):
    if filename.endswith('.bed.gz'):
        input_path = os.path.join(input_folder, filename)
        output_filename = filename.replace('.peaks.bed.gz', '.dedup.peaks.bed.gz')
        output_path = os.path.join(output_folder, output_filename)

        # Read the compressed BED file into a pandas DataFrame
        with gzip.open(input_path, 'rt') as f:
            columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak']
            df = pd.read_csv(f, sep='\t', header=None, names=columns)

        # Group by the first three columns and keep the row with the highest value in the 9th column (qValue)
        df_deduplicated = df.loc[df.groupby(['chrom', 'start', 'end'])['qValue'].idxmax()]

        # Write the deduplicated data to a compressed BED file
        with gzip.open(output_path, 'wt') as f_out:
            df_deduplicated.to_csv(f_out, sep='\t', header=False, index=False)

        print(f"Processed {filename} -> {output_filename}")

print("All files processed and saved in the 'dedup' folder.")

