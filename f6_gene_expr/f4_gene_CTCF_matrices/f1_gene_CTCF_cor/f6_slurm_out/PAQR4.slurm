#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PAQR4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PAQR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAQR5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAQR6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAQR7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAQR8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAQR9 
