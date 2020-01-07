#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CMTM3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CMTM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMTM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMTM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMTM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMTM7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMTM8 
