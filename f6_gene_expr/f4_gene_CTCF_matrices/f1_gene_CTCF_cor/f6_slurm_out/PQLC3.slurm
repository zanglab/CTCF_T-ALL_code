#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PQLC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PQLC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRADC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAM1 
