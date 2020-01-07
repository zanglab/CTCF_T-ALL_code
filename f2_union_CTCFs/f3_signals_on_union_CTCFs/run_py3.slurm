#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=80000
#SBATCH -t 12:00:00
#SBATCH -p gpu 
#SBATCH -A zanglab
#SBATCH -o slurm_py3.out

#Run program
time python py3_QN.py
