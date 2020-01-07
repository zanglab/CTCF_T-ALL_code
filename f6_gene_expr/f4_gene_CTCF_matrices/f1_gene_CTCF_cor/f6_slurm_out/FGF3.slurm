#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FGF3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF8 
