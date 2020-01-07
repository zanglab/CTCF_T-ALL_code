#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CNN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CNN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNNM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNNM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNNM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNNM4 
