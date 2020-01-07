#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ITIH1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ITIH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITIH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITIH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITIH4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITIH5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITIH6 
