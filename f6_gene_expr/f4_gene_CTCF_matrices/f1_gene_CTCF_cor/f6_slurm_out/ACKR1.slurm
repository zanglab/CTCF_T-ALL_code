#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACKR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACKR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACKR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACKR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACKR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACLY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACMSD 
