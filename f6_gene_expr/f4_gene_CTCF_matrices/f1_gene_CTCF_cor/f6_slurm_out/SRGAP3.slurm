#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SRGAP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SRGAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRGN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRMS 
