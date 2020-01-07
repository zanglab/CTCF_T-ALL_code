#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CRLF3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CRLF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRLS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRNKL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRNN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CROCC 
