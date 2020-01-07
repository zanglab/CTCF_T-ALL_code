#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BUD31.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BUD31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BVES 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BYSL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BZW1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BZW2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf105 
