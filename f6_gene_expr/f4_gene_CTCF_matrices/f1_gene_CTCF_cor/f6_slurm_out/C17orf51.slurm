#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C17orf51.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf51 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf64 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf67 
