#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYP3A5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP3A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP3A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP3A7-CYP3A51P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP46A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4A11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4A22 
