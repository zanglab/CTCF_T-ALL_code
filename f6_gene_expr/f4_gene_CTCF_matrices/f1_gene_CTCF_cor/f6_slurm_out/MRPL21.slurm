#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRPL21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL28 
