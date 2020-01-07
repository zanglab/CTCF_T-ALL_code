#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FTH1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FTH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FTH1P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FTL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FTLP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FTMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FTO 
