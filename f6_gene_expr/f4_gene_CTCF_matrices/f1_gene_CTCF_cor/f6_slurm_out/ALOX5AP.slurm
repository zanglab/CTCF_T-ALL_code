#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ALOX5AP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ALOX5AP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALOXE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALPI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALPK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALPK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALPK3 
