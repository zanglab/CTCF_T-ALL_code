#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRTAP19-3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP19-3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP19-8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP2-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP2-3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP2-4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP21-1 
