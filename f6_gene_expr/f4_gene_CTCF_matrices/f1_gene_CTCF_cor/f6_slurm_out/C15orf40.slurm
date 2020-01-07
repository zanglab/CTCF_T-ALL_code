#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C15orf40.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C15orf40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C15orf41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C15orf48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C15orf53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C15orf59 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C15orf61 
