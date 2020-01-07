#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NOTCH2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NOTCH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOTCH2NL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOTCH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOTCH4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOTO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOTUM 
