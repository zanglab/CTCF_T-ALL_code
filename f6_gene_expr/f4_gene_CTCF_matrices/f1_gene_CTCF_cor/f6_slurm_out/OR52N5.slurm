#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR52N5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52N5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52W1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR56A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR56A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR56A4 
