#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLA2G2D.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G2E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G2F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G4B 
