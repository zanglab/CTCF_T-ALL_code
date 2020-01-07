#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SELENBP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SELENBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SELL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SELP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SELPLG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA3B 
