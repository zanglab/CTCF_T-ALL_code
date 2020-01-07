#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR10A3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10AD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10AG1 
