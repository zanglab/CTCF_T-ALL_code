#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C3orf52.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf56 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf67 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf70 
