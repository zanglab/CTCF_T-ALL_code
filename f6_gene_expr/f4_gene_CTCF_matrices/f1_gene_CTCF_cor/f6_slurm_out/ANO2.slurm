#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANO2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANO5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANO6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANO7 
