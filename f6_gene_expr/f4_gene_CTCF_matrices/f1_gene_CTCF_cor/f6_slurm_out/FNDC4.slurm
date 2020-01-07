#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FNDC4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FNDC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FNDC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FNDC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FNDC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FNDC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FNIP1 
