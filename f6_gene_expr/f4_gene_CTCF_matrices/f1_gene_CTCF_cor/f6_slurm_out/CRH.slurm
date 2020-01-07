#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CRH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CRH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRHBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRHR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRHR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRIM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRIP1 
