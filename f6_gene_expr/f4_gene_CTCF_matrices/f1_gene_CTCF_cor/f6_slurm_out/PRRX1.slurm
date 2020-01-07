#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRRX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRRX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRRX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS2 
