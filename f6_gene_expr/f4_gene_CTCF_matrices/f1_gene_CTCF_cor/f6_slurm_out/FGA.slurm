#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FGA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FGA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGD4 
