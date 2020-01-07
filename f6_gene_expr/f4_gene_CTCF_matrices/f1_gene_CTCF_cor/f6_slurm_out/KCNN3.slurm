#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNN3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNQ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNQ2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNQ3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNQ4 
