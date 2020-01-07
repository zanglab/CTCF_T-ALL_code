#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MTURN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MTURN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTUS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTUS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTX3 
