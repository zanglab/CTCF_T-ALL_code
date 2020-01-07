#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o INTS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i INTS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INTS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INTS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INTS4P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INTS5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INTS6 
