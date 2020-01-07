#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABCF3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCG8 
