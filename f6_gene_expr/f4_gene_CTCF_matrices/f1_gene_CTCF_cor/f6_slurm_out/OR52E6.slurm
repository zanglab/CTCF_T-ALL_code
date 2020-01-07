#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR52E6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52E6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52E8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52I1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52I2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52K1 
