#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMCHD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMCHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMCO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMCO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMCO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMCO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMCP 
