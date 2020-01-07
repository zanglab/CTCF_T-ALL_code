#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SAMM50.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SAMM50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAMSN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAP130 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAP18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAP25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAP30 
