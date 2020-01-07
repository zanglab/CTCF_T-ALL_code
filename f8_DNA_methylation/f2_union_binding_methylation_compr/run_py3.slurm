#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p gpu 
#SBATCH -A zanglab 
#SBATCH -o run3.out

time python py3_methy_compr_add_UnCoScored_regions.py
