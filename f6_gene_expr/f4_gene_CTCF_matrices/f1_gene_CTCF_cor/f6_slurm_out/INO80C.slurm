#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o INO80C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i INO80C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INO80D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INO80E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INPP4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INPP4B 
