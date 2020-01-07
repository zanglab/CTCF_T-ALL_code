#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PSD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PSD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSEN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSEN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSENEN 
