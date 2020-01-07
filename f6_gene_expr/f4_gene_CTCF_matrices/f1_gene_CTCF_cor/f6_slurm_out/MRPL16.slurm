#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRPL16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL20 
