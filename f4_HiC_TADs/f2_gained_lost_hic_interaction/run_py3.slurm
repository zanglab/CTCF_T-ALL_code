#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p gpu 
#SBATCH -A zanglab
#SBATCH -o slurm_py3.out
time python 3_compr_paired_ttest.py
