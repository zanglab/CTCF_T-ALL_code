#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNMB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNMB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNMB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNMB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNMB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNN2 
