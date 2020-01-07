#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NCKAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NCKAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCKAP1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCKAP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCKAP5L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCKIPSD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCL 
