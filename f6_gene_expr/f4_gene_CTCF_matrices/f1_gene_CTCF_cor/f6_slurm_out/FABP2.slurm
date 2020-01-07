#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FABP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FABP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FABP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FABP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FABP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FABP5P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FABP6 
