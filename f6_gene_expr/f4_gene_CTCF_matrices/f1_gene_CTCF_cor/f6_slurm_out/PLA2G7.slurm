#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLA2G7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2G7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLA2R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLAA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLAC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLAC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLAC8L1 
