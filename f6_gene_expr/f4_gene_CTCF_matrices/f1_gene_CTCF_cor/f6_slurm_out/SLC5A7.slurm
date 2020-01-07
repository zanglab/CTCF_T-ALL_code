#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC5A7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC5A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC5A8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC5A9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC6A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC6A10P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC6A11 
