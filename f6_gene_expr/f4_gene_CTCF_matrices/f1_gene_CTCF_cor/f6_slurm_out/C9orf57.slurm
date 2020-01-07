#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C9orf57.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf57 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf64 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf66 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf72 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf78 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf84 
