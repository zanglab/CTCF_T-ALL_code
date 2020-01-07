#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPS4Y1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS4Y1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS4Y2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS6KA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS6KA2 
