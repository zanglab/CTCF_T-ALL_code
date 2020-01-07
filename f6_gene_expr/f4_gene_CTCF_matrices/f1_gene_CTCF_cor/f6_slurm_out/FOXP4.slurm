#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FOXP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXQ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXRED1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXRED2 
