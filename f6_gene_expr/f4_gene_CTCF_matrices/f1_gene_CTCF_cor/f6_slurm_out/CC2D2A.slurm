#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CC2D2A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CC2D2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CC2D2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCAR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCBE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC102A 
