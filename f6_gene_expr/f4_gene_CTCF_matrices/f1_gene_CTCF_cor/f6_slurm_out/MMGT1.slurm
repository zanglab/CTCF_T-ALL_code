#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MMGT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MMGT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP13 
