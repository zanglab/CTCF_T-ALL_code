#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SIVA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SIVA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIX5 
