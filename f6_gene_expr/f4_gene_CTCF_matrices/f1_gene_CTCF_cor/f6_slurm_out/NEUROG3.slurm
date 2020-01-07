#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NEUROG3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NEUROG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEXN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NF1P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFAM1 
