#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLEK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHA5 
