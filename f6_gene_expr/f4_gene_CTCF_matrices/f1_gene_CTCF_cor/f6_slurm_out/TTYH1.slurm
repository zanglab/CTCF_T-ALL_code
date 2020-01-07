#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TTYH1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TTYH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTYH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTYH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBA1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBA1B 
