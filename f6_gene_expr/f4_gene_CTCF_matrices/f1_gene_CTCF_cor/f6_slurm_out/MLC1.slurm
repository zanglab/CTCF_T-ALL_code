#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MLC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MLC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLEC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLH3 
