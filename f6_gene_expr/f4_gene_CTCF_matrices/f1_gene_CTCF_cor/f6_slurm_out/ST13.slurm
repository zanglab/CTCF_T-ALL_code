#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ST13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ST13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST13P4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST20-MTHFS 
