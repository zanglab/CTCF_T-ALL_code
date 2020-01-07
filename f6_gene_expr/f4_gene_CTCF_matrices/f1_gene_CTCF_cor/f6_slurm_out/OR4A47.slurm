#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR4A47.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4A47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4C6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4D10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4D5 
