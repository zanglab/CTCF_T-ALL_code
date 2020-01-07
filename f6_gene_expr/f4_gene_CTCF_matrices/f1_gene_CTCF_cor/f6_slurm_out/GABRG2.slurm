#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GABRG2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRQ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRR2 
