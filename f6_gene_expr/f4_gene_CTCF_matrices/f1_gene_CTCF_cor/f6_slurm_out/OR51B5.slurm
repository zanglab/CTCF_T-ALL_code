#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR51B5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51B5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51B6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51E1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51E2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51F1 
