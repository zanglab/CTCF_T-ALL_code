#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C12orf29.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf45 
