#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRTAP4-9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP4-9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP5-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP5-10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP5-11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP5-2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP5-3 
