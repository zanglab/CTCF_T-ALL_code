#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TEX261.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX261 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX264 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX33 
