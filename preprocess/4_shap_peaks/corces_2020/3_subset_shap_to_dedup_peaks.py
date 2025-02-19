import os
import gzip
import numpy as np
import deepdish
import shutil
import pandas as pd
from multiprocessing import Pool


def subset_mean_files(args):

    cluster, input_dir, dedup_peaks_dir, output_dir, shap_type = args

    """
    Subset the mean .h5 and .npz files using dedup peaks and save to the output directory.
    """
    print(f"Processing cluster: {cluster}")

    # Define file paths
    mean_h5_file = os.path.join(input_dir, cluster, "mean", f"{cluster}.mean.original_peaks.{shap_type}_scores.{shap_type}_scores.h5")
    mean_npz_file = os.path.join(input_dir, cluster, "mean", f"{cluster}.mean.original_peaks.{shap_type}_scores.{shap_type}_scores.npz")
    mean_bed_file = os.path.join(input_dir, cluster, "mean", f"{cluster}.mean.original_peaks.{shap_type}_scores.interpreted_regions.bed")
    dedup_peaks_file = os.path.join(dedup_peaks_dir, f"{cluster}.overlap.dedup.peaks.bed.gz")
    output_h5_file = os.path.join(output_dir, cluster, "mean", f"{cluster}.mean.specific_peaks.{shap_type}_scores.{shap_type}_scores.h5")
    output_npz_file = os.path.join(output_dir, cluster, "mean", f"{cluster}.mean.specific_peaks.{shap_type}_scores.{shap_type}_scores.npz")
    output_bed_file = os.path.join(output_dir, cluster, "mean", f"{cluster}.mean.specific_peaks.{shap_type}_scores.interpreted_regions.bed")

    # Create output directories if not exists
    os.makedirs(os.path.dirname(output_h5_file), exist_ok=True)

    dedup_peaks = pd.read_table(dedup_peaks_file, header=None, usecols=[0, 1, 2, 3], dtype=str)
    mean_bed = pd.read_table(mean_bed_file, header=None, dtype=str)

    # Find matching rows and their indices
    mean_bed['index'] = mean_bed.index
    merged = pd.merge(mean_bed, dedup_peaks, on=[0, 1, 2, 3], how='inner')
    matching_indices = merged['index'].tolist()

    if not matching_indices:
        print(f"No matching peaks found for cluster {cluster}. Skipping.")
        return

    # Load the mean .h5 and .npz files
    data_h5 = deepdish.io.load(mean_h5_file)

    # Subset the SHAP scores using matching indices
    subset_shap = data_h5['shap']['seq'][matching_indices]
    subset_projected_shap = data_h5['projected_shap']['seq'][matching_indices]
    subset_onehot = data_h5['raw']['seq'][matching_indices]

    # Save subset data to output files
    deepdish.io.save(output_h5_file, {
        'raw': {'seq': subset_onehot},
        'shap': {'seq': subset_shap},
        'projected_shap': {'seq': subset_projected_shap}
    }, compression='blosc')

    np.savez(output_npz_file, subset_shap)

    # Save the subsetted bed file
    merged.to_csv(output_bed_file, sep='\t', header=False, index=False)
    print(f"Cluster {cluster} processed and saved.")


def process_clusters(input_dir, dedup_peaks_dir, output_dir, shap_types, max_processes=20):
    """
    Iterate through clusters and process the mean files.
    """

    args_list = []

    for cluster in os.listdir(input_dir):
        for shap_type in shap_types:
            args_list.append((cluster, input_dir, dedup_peaks_dir, output_dir, shap_type))

    # Use Pool to parallelize the processing
    with Pool(processes=max_processes) as pool:
        pool.map(subset_mean_files, args_list)


# Directories and parameters
input_dir = "/oak/stanford/groups/akundaje/projects/neuro-variants/peak_shap/corces_2020/original_peaks"
dedup_peaks_dir = "/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/corces_2020/peaks/dedup_overlap"
output_dir = "/oak/stanford/groups/akundaje/projects/neuro-variants/peak_shap/corces_2020/specific_peaks"
shap_types = ['counts']  # Add 'profile' if needed

# Run the processing
process_clusters(input_dir, dedup_peaks_dir, output_dir, shap_types, max_processes=10)
