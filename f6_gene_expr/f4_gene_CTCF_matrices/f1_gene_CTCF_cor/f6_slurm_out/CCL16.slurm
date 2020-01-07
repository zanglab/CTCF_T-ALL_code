#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCL16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL20 
