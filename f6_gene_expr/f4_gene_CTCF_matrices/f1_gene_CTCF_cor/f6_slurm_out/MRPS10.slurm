#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRPS10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS16 
