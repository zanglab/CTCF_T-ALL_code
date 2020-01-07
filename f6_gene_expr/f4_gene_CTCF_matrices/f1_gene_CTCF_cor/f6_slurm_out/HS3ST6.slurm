#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HS3ST6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HS3ST6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HS6ST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HS6ST2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HS6ST3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSBP1L1 
