#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TUBA8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBAL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBB2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBB2B 
