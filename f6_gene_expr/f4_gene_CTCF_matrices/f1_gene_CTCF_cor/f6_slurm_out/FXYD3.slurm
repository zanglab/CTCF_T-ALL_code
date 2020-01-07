#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FXYD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FXYD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FXYD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FXYD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FXYD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FXYD6-FXYD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FXYD7 
