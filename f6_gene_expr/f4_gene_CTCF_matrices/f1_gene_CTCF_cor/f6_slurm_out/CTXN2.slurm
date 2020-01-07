#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CTXN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CTXN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTXN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUBN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUEDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUEDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUL1 
