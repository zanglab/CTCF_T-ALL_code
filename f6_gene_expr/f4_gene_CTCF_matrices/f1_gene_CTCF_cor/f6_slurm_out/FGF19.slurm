#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FGF19.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF23 
