#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o py1.out

#time python 1_cancertype_specific_binding_ATAC_sig_abstract.py

time python 1_cancertype_specific_binding_ATAC_sig_abstract.py  > py1.out
time sh py1.out
