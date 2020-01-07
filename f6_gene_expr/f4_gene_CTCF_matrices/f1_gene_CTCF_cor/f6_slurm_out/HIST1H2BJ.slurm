#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIST1H2BJ.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BJ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BO 
