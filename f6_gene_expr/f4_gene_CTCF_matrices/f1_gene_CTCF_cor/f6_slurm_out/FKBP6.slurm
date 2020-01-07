#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FKBP6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP9P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBPL 
