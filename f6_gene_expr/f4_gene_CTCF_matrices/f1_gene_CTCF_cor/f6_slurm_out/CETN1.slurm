#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CETN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CETN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CETN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CETN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CETN4P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CETP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP126 
