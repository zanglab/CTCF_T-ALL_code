#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ASB6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB9P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASCC1 
