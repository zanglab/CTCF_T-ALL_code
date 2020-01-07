#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCL21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL26 
