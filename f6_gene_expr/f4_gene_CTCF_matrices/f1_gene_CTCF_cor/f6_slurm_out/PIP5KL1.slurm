#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIP5KL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIP5KL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIPOX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIPSL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIRT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PISD 
