#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABHD4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABI2 
