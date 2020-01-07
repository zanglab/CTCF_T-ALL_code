#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TEX35.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX9 
