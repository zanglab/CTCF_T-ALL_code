#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TGM3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TGM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGM7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGOLN2 
