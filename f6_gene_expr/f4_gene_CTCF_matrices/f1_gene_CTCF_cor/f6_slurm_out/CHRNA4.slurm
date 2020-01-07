#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHRNA4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNA9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNB1 
