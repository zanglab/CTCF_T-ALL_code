#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POU4F2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POU4F2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU4F3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU5F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU5F1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU5F1P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU5F1P4 
