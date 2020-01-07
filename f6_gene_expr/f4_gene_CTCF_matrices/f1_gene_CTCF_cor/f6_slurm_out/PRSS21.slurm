#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRSS21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS30P 
