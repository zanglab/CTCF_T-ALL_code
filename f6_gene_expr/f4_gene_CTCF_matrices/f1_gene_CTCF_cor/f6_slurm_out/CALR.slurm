#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CALR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CALR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CALR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CALU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CALY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMK1D 
