#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SGSM2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SGSM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGSM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGTA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGTB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH2B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH2B2 
