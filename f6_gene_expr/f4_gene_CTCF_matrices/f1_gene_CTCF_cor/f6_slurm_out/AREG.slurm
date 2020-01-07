#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AREG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AREG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AREL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARF5 
