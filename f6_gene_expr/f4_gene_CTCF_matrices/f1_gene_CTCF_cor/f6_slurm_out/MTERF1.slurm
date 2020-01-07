#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MTERF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MTERF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTERF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTERF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTERF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTF2 
