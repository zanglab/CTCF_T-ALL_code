#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HNF1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HNF1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNF1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNF4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNF4G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPA0 
