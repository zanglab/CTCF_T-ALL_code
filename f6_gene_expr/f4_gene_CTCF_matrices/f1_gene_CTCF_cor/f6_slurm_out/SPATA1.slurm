#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPATA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA18 
