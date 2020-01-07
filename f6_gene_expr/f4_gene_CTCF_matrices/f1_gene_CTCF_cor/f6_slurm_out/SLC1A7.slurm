#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC1A7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC1A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC20A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC20A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A11 
