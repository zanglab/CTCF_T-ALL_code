#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SERINC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SERINC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERINC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERINC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERINC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERINC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERP1 
