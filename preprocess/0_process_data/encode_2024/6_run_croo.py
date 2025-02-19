import os
import json
import subprocess

metadata_json_dir = 'pipeline_metadata_jsons'
croo_dir = 'gs://neuro-variants/data/processed/encode_2024/croo'

for file_name in os.listdir(metadata_json_dir):
    sample_name = '.'.join(file_name.split('.')[:4])
    print(sample_name)

    cmd = 'croo ' + os.path.join(metadata_json_dir, file_name) \
          + ' --out-dir ' + os.path.join(croo_dir, sample_name) \
          + ' --method copy'
    cmd = cmd.split()
    print(cmd)
    ret = subprocess.call(cmd)

