#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPINK9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINK9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPIRE1 
