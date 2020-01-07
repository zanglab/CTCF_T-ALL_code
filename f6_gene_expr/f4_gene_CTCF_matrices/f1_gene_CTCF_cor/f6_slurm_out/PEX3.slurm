#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PEX3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX5L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PF4 
