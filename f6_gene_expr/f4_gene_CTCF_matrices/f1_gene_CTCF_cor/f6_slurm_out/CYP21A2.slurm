#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYP21A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP21A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP24A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP26A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP26B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP26C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP27A1 
