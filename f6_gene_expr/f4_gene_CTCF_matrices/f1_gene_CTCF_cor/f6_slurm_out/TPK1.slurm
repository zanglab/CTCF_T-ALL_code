#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TPK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TPK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPM3P9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPM4 
