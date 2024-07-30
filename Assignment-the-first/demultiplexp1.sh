#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-user=wesg@uoregon.edu
#SBATCH --mail-type=END

conda activate bgmp_py.mplib
conda list

/usr/bin/time -v ./mean_qual_fq.py 

conda deactivate
