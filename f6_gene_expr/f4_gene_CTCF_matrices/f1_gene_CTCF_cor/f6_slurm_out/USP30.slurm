#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o USP30.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i USP30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP32P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP32P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP33 
