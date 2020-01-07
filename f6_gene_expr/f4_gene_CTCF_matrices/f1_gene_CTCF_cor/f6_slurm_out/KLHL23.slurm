#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLHL23.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL29 
