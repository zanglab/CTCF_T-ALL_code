#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM222B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM222B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM227A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM227B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM228A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM228B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM229A 
