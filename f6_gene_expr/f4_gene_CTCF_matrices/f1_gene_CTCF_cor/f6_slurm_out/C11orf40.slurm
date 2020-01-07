#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C11orf40.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf49 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf53 
