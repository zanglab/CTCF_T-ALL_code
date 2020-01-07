#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC9A4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9A7P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9A8 
