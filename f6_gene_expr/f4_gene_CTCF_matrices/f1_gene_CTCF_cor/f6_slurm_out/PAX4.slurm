#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PAX4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PAX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAX7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAX8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAX9 
