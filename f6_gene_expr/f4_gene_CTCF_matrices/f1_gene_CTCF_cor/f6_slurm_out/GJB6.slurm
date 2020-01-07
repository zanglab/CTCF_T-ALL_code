#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GJB6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GJB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJB7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJD2 
