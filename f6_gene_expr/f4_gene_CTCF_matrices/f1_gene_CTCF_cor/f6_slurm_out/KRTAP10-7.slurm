#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRTAP10-7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP10-7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP10-8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP10-9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP11-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP12-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP12-2 
