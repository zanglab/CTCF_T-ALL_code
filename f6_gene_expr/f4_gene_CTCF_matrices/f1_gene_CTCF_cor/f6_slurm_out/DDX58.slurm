#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DDX58.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX59 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX60 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX60L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEAF1 
