#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABCA11P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA11P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA17P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA3 
