#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRTAP5-4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP5-4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP5-5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP5-6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP5-7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP5-8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP5-9 
