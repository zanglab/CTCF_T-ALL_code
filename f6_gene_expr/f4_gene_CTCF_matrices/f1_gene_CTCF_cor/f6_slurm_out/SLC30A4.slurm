#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC30A4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC30A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC30A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC30A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC30A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC30A8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC30A9 
