#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FASN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FASN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FASTK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FASTKD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FASTKD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FASTKD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FASTKD5 
