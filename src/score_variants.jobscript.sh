#!/bin/bash

set -e
set -u
set -o pipefail

module load cuda/11.2
module load cudnn/8.1

# https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
eval "$(conda shell.bash hook)"
conda activate chrombpnet

echo "Live"
python "$@"

