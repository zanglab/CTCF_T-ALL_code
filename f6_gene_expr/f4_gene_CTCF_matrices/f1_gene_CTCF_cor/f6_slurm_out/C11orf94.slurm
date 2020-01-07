#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C11orf94.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf94 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf95 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf96 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf97 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf98 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf10 
