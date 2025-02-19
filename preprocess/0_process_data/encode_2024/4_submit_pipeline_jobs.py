import os
import subprocess

input_json_dir = 'pipeline_input_jsons'
wdl_file = '/opt/atac-seq-pipeline/atac.wdl'
output_dir = 'gs://neuro-variants/data/processed/atac_pipeline_outputs'
port = '8090'

for file_name in os.listdir(input_json_dir):
    sample_name = '.'.join(file_name.split('.')[:4])
    print(sample_name)

    cmd = 'caper submit ' + wdl_file + ' -i ' + input_json_dir \
          + '/' + file_name + ' -s ' + sample_name \
          + ' --gcp-loc-dir=' + output_dir + '/.caper_tmp --port ' \
          + port
    cmd = cmd.split()
    print(cmd)
    ret = subprocess.call(cmd)

