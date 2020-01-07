#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AURKC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AURKC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AUTS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AVEN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AVIL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AVL9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AVP 
