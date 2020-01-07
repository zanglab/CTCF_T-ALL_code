#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRPL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL35 
