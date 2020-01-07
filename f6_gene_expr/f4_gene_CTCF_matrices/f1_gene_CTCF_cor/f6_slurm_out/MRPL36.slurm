#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRPL36.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL40 
