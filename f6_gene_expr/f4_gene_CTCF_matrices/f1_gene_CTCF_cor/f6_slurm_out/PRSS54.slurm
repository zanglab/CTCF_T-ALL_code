#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRSS54.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS55 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS56 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS57 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS8 
