#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PUS10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PUS10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PUS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PUS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PUS7L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PUSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PVALB 
