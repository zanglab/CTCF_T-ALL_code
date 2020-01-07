#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WNT11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT3A 
