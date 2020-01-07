#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DFFB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DFFB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGAT2L6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGCR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGCR5 
