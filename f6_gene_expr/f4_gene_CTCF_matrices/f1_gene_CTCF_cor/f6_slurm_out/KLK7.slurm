#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLK7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLK7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLK8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLK9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLKB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLKP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLLN 
