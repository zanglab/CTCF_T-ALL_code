#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRTAP12-3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP12-3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP13-2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP13-4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP16-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP17-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP19-1 
