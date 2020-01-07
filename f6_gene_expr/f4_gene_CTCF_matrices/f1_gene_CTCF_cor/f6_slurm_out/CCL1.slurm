#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL15-CCL14 
