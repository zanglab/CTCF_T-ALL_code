#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PHKG2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PHKG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHLDA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHLDA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHLDA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHLDB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHLDB2 
