#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAP4K3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP4K3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP4K4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP4K5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP6D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP7 
