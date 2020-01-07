#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABCB7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCB7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCB8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCB9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCC10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCC11 
