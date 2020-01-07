#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYO7A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO9A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO9B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYOC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYOCD 
