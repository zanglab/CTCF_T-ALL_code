#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TP53BP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53BP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53I11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53I13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53I3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53INP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53INP2 
