#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CFAP99.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP99 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFC1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFDP1 
