#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CMKLR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CMKLR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMPK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMPK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMSS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMTM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMTM2 
