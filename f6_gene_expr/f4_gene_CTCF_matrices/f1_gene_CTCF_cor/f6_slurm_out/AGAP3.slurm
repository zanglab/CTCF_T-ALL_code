#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AGAP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AGAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGAP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGAP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGAP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGAP7P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGAP9 
