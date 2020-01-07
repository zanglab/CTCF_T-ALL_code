#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRTAP10-12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP10-12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP10-2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP10-3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP10-4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP10-5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP10-6 
