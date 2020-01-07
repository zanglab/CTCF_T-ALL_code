#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABCD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCF2 
