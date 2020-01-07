#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DDX43.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX49 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX50 
