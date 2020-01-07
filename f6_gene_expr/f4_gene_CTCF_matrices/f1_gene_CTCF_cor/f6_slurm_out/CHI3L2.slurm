#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHI3L2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHI3L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHIA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHIAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHIC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHIC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHID1 
