#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EVA1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EVA1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EVA1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EVC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EVC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EVI2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EVI2B 
