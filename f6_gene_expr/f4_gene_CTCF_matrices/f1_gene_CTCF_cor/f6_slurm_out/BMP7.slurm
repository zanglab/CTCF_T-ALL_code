#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BMP7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BMP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMP8A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMP8B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMPER 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMPR1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMPR1B 
