#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C6orf48.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf89 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C7 
