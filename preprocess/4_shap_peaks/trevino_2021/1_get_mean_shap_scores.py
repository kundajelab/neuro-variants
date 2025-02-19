import h5py
import numpy as np
import os
import deepdish
import shutil
from multiprocessing import Pool


def get_mean_shap(args):
    peak_type, cluster, shap_type, shap_dir = args
    print(f'Processing: {peak_type}, {cluster}, {shap_type}')

    outdir = os.path.join(shap_dir, peak_type, cluster, 'mean')
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    outfile_npz = os.path.join(outdir, '.'.join([cluster, 'mean', peak_type, shap_type + '_scores',
                                                 shap_type + '_scores', 'npz']))
    outfile_h5 = os.path.join(outdir, '.'.join([cluster, 'mean', peak_type, shap_type + '_scores',
                                                shap_type + '_scores', 'h5']))

    # if not os.path.isfile(outfile_npz):
    shap_dict = {}
    projected_shap_dict = {}
    onehot = []

    # Processing folds
    for fld in range(5):
        fold = 'fold_' + str(fld)
        print(fold)

        infile = os.path.join(shap_dir, peak_type, cluster, fold,
                                '.'.join([cluster, fold, peak_type, shap_type + '_scores',
                                        shap_type + '_scores', 'h5']))
        shap_dict[fold] = deepdish.io.load(infile, '/shap/seq')
        projected_shap_dict[fold] = deepdish.io.load(infile, '/projected_shap/seq')
        if fold == 'fold_0':
            onehot = deepdish.io.load(infile, '/raw/seq')

    mean_shap = np.mean(np.array([shap_dict[fold] for fold in shap_dict]), axis=0)
    mean_projected_shap = np.mean(np.array([projected_shap_dict[fold] for fold in projected_shap_dict]), axis=0)
    d = {
        'raw': {'seq': onehot},
        'shap': {'seq': mean_shap},
        'projected_shap': {'seq': mean_projected_shap}
    }

    np.savez(outfile_npz, mean_shap)
    deepdish.io.save(outfile_h5, d, compression='blosc')

    onehot_file = os.path.join(outdir, '.'.join([cluster, 'mean', peak_type, shap_type + '_scores',
                                                 'onehot', 'npz']))
    np.savez(onehot_file, onehot)

    shutil.copy(os.path.join(shap_dir, peak_type, cluster, 'fold_0',
                                '.'.join([cluster, 'fold_0', peak_type, shap_type + '_scores',
                                        'interpreted_regions', 'bed'])),
                os.path.join(outdir,
                                '.'.join([cluster, 'mean', peak_type, shap_type + '_scores',
                                        'interpreted_regions', 'bed'])))

def main_parallel_processing(shap_dir, peak_types, shap_types, max_processes=20):
    # Create a list of arguments to be passed to the worker function
    args_list = []
    for peak_type in peak_types:
        print(peak_type)
        print()

        # for cluster in os.listdir(os.path.join(shap_dir, peak_type)):
        for cluster in ['trevino_2021.c14']:
            print(cluster)
            print()

            for shap_type in shap_types:
                args_list.append((peak_type, cluster, shap_type, shap_dir))

    # Use Pool to parallelize the processing
    with Pool(processes=max_processes) as pool:
        pool.map(get_mean_shap, args_list)


shap_dir = '/oak/stanford/groups/akundaje/projects/neuro-variants/peak_shap/trevino_2021'
shap_types = ['counts'] #, 'profile']
peak_types = ['original_peaks'] #, 'specific_peaks']

main_parallel_processing(shap_dir, peak_types, shap_types, max_processes=10)

