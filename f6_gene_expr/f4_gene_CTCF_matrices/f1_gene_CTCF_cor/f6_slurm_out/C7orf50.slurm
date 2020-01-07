#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C7orf50.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C7orf50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C7orf55-LUC7L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C7orf57 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C7orf61 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C7orf77 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8A 
