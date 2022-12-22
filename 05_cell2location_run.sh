#!/usr/bin/bash
#SBATCH --job-name 1-1
#SBATCH --account PCON0022
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=50GB
#SBATCH --gpus-per-node=1

module load cuda
cd /fs/ess/PCON0022/guoqi/NC-snrna/Spatial/Deconvolution_6scRNA

python cell2location_pipeline.py --sample 1-1

python cell2location_pipeline.py --sample 1-7

python cell2location_pipeline.py --sample 18-64

python cell2location_pipeline.py --sample 2-10

python cell2location_pipeline.py --sample 2-3

python cell2location_pipeline.py --sample 2-5

python cell2location_pipeline.py --sample 2-8

python cell2location_pipeline.py --sample T4857