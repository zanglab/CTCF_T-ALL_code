#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC1A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC1A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC1A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC1A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC1A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC1A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC1A6 
