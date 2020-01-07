#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NPTXR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NPTXR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPW 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPY1R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPY2R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPY4R 
