#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SYNE1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNE4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNGAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNGR1 
