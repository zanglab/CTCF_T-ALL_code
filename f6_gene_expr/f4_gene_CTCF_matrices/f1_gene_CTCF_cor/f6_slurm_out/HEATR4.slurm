#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HEATR4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HEATR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEATR5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEATR5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEATR6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEATR9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEBP1 
