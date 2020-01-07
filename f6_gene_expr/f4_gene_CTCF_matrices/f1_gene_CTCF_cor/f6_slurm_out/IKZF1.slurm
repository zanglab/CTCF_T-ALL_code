#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IKZF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IKZF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IKZF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IKZF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IKZF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IKZF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL10 
