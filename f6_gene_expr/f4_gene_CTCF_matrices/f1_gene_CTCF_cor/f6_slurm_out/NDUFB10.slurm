#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NDUFB10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFB10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFB11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFB5 
