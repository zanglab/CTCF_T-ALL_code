#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC6A5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC6A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC6A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC6A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC6A8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC6A9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A1 
