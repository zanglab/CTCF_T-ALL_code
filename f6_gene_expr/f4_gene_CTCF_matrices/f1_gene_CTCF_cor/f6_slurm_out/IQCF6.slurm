#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IQCF6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCJ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCJ-SCHIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCK 
