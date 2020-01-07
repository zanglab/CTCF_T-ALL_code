#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=80000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab

time python py3_union_all_peak_summits.py
