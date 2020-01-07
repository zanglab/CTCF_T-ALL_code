#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EBF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EBF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EBF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EBF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EBF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EBI3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EBLN1 
