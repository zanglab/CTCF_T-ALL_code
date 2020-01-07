#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PKP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PKP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G12A 
