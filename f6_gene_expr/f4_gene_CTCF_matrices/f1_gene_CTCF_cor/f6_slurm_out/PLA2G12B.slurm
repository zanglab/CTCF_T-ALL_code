#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLA2G12B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G12B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G2C 
