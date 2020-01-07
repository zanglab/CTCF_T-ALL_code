#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLA2G4C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G4C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G4D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G4E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G4F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G6 
