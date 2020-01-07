#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p gpu
#SBATCH -A zanglab

#Run program
time python 1_TALL_gained_lost_binding_signal_compr.py

