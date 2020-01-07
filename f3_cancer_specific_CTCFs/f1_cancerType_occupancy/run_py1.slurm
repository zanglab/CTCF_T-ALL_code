#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p gpu 
#SBATCH -A zanglab
#SBATCH -o slurm_py1.out
#Run program
time python 1_binding_occupancy_sep_cancerType.py
