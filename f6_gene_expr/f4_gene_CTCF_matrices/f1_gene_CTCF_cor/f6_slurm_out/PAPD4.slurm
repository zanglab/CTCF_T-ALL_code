#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PAPD4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PAPD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAPD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAPD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAPLN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAPOLA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAPOLB 
