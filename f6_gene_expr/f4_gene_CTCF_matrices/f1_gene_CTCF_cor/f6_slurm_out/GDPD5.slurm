#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GDPD5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GDPD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDPGP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GEM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GEMIN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GEMIN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GEMIN5 
