#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KDM6B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KDM6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDM7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDM8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDSR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KEAP1 
