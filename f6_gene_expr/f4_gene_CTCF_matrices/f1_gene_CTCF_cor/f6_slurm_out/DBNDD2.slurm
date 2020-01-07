#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DBNDD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DBNDD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DBNL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DBR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DBT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DBX1 
