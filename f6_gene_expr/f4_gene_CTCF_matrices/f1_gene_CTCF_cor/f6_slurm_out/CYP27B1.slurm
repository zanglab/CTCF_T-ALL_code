#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYP27B1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP27B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP27C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2A13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2B6 
