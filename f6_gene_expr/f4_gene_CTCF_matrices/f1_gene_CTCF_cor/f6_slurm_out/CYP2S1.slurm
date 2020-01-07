#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYP2S1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2S1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2U1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2W1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP39A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP3A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP3A43 
