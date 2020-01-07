#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAP3K7CL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K7CL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP4K1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP4K2 
