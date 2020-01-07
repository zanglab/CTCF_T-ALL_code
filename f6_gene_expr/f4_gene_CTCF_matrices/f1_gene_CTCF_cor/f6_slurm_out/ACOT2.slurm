#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACOT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOT7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOT8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOT9 
