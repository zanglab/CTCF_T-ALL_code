#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYP2B7P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2B7P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2C18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2C19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2C8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2C9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2D6 
