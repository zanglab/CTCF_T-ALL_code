#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NPL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPLOC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPNT 
