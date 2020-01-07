#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o S100A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i S100A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i S100A10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i S100A11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i S100A12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i S100A13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i S100A14 
