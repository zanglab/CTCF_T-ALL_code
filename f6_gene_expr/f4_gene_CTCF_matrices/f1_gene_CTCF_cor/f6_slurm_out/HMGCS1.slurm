#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HMGCS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGCS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGCS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGN2P46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGN3 
