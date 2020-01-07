#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o slurm_py1.out3

#Run program
time python 1_ttest_on_signal_sep_cancerType.py -i CRC