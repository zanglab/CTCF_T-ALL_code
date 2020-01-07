#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIP4K2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIP4K2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIP4K2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIP5K1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIP5K1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIP5K1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIP5K1P1 
