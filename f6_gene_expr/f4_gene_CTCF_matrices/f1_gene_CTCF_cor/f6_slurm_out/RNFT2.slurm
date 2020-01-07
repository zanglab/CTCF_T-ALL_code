#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNFT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNFT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNGTT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNLS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNPC3 
