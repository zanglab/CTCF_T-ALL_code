#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FOXD4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXD4L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXD4L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXD4L4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXD4L5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXD4L6 
