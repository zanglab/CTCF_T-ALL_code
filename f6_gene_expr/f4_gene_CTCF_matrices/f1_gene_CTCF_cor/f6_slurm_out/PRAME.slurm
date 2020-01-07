#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRAME.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAME 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAMEF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAMEF12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAMEF13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAMEF14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAMEF15 
