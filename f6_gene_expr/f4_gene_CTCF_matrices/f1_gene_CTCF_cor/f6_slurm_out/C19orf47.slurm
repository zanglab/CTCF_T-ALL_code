#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C19orf47.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C19orf47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C19orf48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C19orf53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C19orf54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C19orf57 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C19orf66 
