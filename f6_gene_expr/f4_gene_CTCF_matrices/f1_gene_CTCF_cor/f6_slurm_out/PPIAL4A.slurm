#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPIAL4A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIAL4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIAL4C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIAL4D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIAL4E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIAL4F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIAL4G 
