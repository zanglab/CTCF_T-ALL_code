#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C17orf75.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf75 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf78 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf80 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf97 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf98 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf99 
