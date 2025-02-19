import h5py
import numpy as np
import os
import deepdish
import shutil
from multiprocessing import Pool


def get_mean_shap(args):
    cluster, shap_type, shap_dir = args
    print(f'Processing: {cluster}, {shap_type}')

    outdir = os.path.join(shap_dir, cluster, 'mean')
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    outfile_h5 = os.path.join(outdir, '.'.join([cluster, 'mean', 'rare', 'shap', 'variant_shap',
                                                shap_type, 'h5']))

    # if not os.path.isfile(outfile_npz):
    shap_dict = {}
    projected_shap_dict = {}
    onehot = []
    variant_ids = []
    alleles = []

    # Processing folds
    for fld in range(5):
        fold = 'fold_' + str(fld)
        print(fold)

        infile = os.path.join(shap_dir, cluster, fold,
                                '.'.join([cluster, fold, 'rare', 'shap', 'variant_shap',
                                        shap_type, 'h5']))
        shap_dict[fold] = deepdish.io.load(infile, '/shap/seq')
        projected_shap_dict[fold] = deepdish.io.load(infile, '/projected_shap/seq')
        if fold == 'fold_0':
            onehot = deepdish.io.load(infile, '/raw/seq')
            variant_ids = deepdish.io.load(infile, '/variant_ids')
            alleles = deepdish.io.load(infile, '/alleles')

    mean_shap = np.mean(np.array([shap_dict[fold] for fold in shap_dict]), axis=0)
    mean_projected_shap = np.mean(np.array([projected_shap_dict[fold] for fold in projected_shap_dict]), axis=0)
    d = {
        'raw': {'seq': onehot},
        'shap': {'seq': mean_shap},
        'projected_shap': {'seq': mean_projected_shap},
        'variant_ids': variant_ids,
        'alleles': alleles
    }

    deepdish.io.save(outfile_h5, d, compression='blosc')


def main_parallel_processing(shap_dir, shap_types, max_processes=20):
    # Create a list of arguments to be passed to the worker function
    args_list = []
    for cluster in os.listdir(shap_dir):
        print(cluster)
        print()

        for shap_type in shap_types:
            args_list.append((cluster, shap_type, shap_dir))

    # Use Pool to parallelize the processing
    with Pool(processes=max_processes) as pool:
        pool.map(get_mean_shap, args_list)


shap_dir = '/oak/stanford/groups/akundaje/projects/neuro-variants/variant_shap/rare/K562_bias/trevino_2021'
shap_types = ['counts'] #, 'profile']

main_parallel_processing(shap_dir, shap_types, max_processes=10)

