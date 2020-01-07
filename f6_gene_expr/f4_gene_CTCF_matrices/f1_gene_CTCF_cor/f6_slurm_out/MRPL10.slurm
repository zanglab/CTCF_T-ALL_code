#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRPL10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL15 
