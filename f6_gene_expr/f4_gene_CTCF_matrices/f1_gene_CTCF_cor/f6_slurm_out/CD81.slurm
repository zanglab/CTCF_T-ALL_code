#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD81.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD81 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD82 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD83 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD84 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD86 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD8A 
