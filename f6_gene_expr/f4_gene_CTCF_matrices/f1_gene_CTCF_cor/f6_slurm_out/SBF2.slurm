#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SBF2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SBF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SBK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SBK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SBK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SBNO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SBNO2 
