#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIN1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIN4 
