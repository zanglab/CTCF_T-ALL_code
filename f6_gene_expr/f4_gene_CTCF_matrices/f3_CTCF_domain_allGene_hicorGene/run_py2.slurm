#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o out_py2.log
#Run program
time python 2_CTCF_domain_NoLimit_allGene_hicorGene.py
