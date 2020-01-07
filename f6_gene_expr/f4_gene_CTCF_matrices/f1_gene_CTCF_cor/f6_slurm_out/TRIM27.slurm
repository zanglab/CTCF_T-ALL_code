#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRIM27.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM32 
