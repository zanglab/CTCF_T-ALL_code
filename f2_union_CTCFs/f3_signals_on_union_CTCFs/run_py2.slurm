#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=60000
#SBATCH -t 6:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o slurm_py2.out
#Run program

time python py2_combine_all_csv.py
