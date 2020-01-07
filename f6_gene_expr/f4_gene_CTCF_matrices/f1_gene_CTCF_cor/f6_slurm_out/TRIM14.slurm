#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRIM14.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM16L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM2 
