#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNX5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOAT1 
