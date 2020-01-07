#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATXN10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATXN10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATXN1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATXN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATXN2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATXN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATXN3L 
