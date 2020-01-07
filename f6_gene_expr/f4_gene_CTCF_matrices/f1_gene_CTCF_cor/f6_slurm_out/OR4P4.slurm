#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR4P4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4P4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4Q3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51B4 
