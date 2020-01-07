#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CNGA3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CNGA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNGA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNGB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNGB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNIH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNIH2 
