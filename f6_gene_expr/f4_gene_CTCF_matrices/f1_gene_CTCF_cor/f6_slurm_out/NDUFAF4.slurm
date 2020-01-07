#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NDUFAF4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFAF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFAF4P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFAF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFAF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFAF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFB1 
