#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C11orf74.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf74 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf80 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf86 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf87 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf88 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf91 
