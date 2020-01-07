#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLOD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLOD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLRG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLS3 
