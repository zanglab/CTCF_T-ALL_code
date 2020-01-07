#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRTFDC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRTFDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRTG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRTN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRUNE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRY 
