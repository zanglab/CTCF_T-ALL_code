#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STXBP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STXBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STXBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STXBP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STXBP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STXBP5L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STXBP6 
