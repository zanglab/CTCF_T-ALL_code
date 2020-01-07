#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRPL41.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL45P2 
