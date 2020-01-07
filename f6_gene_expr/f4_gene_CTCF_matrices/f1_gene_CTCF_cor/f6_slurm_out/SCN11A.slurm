#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SCN11A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SCN11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCN1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCN1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCN2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCN2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCN3A 
