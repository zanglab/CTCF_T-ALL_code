#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WBP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WBP2NL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WBP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDFY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDFY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDFY3 
