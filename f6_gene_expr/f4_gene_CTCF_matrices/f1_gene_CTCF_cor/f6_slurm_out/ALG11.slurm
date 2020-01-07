#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ALG11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG1L2 
