#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WNT4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT7B 
