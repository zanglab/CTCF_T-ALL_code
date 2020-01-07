#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HARS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HARS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAUS1 
