#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRTAP4-1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP4-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP4-11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP4-12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP4-4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP4-6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP4-8 
