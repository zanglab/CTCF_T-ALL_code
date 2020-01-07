#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCDH11Y.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDH11Y 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDH12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDH15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDH17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDH18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDH19 
