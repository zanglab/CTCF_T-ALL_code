#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 4:00:00
#SBATCH -p gpu
#SBATCH -A zanglab

#Run program
time python 3_cancerType_gained_constitutive_signal_compr.py
