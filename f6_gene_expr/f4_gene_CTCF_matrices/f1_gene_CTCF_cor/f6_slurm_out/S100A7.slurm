#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o S100A7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i S100A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i S100A7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i S100A8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i S100A9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i S100B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i S100G 
