#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD8B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD8B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD93 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD96 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD99 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD99L2 
