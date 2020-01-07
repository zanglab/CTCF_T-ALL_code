#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DDX39B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX39B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX3X 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX3Y 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX42 
