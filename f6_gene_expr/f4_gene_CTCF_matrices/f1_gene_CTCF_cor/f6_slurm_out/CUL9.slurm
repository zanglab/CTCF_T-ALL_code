#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CUL9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CUL9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUTA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUTC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUZD1 
