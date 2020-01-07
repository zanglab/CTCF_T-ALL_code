#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDKN2D.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKN2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDNF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDON 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDPF1 
