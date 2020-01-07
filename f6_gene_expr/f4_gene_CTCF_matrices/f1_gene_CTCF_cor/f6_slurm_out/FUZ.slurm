#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FUZ.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FUZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FXN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FXR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FXR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FXYD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FXYD2 
