#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GJA9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GJA9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJB5 
