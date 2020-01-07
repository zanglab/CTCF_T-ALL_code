#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TIFA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TIFA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIFAB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIGAR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIGD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIGD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIGD3 
