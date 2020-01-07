#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o py2.out

#time python 1_cancertype_specific_binding_ATAC_sig_abstract.py
time python 2_cancertype_specific_binding_ATAC_sig_abstract_panCancer_saved.py > py2.out
time sh py2.out
