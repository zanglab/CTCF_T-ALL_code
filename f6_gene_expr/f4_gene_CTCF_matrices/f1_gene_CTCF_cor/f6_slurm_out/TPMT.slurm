#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TPMT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TPMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPPP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPPP2 
