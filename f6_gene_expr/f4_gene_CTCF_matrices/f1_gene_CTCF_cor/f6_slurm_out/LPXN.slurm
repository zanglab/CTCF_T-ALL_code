#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LPXN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LPXN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRBA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRCH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRCH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRCH3 
