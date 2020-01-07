#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TUBB3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBB4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBB4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBB8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBBP5 
