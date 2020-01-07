#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab

#Run program
time python py1_extract_gained_lost_binding_pattern_rerank.py
