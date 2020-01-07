#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CAPN6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPNS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPNS2 
