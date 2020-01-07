#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ASB17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB5 
