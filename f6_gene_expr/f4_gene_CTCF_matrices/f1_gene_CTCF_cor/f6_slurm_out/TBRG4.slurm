#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TBRG4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TBRG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBX10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBX15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBX18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBX19 
