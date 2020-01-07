#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CPSF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CPSF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPSF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPSF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPSF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPSF4L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPSF6 
