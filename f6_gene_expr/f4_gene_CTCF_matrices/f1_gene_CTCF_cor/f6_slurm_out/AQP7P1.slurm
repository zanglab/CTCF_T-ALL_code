#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AQP7P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP7P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP7P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AR 
