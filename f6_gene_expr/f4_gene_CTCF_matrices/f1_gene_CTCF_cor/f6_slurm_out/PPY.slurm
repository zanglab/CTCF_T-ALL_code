#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPY.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPY2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PQBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PQLC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PQLC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PQLC2L 
