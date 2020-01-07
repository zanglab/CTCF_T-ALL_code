#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ODAM.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ODAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ODC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ODF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ODF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ODF2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ODF3 
