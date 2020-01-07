#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRIP12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIP12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIP13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIQK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIT1 
