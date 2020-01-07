#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIAA1109.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1109 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1143 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1147 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1191 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1210 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1211 
