#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GRHPR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GRHPR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRID1 
