#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NDUFB6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFB7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFB8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFB9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFC2 
