#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MFHAS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MFHAS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFNG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFRP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFSD1 
