#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYP4B1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4F11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4F12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4F2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4F22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4F24P 
