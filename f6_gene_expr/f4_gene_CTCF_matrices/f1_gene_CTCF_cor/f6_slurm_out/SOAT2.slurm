#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SOAT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SOAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOCS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOCS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOCS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOCS4 
