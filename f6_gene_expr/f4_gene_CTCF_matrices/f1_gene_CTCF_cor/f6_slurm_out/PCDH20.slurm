#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCDH20.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDH20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDH7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDH8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDH9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHA10 
