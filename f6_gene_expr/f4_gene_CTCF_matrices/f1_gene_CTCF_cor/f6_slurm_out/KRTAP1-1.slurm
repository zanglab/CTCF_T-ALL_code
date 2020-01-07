#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRTAP1-1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP1-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP1-3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP1-5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP10-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP10-10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP10-11 
