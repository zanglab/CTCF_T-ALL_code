#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o py2.out
#Run program
time python 2_T-ALL_gained_lost_notch_hic_interaction.py
