#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GABRB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRG1 
