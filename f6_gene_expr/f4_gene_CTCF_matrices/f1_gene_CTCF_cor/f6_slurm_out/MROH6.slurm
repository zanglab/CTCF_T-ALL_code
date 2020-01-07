#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MROH6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MROH6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MROH7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MROH7-TTC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MROH8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MROH9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL1 
