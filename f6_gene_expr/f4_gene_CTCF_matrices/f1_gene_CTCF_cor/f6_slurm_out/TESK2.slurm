#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TESK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TESK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TESPA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TET1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TET2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TET3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX10 
