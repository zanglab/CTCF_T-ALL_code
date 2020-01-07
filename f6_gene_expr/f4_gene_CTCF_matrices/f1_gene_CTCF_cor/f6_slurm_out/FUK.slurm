#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FUK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FUK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUNDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUNDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUNDC2P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUOM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FURIN 
