#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRRF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRRF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRTO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRVI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A10 
