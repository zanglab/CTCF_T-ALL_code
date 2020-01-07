#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o out_py1.log
#Run program
time python 1_CTCF_domain_GT100K_LT1M_allGene_hicorGene.py
