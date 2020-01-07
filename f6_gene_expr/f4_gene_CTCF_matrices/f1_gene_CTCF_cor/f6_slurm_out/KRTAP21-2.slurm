#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRTAP21-2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP21-2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP26-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP29-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP3-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP3-2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP3-3 
