#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPRYD7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPRYD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPSB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPSB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPSB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPSB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPTA1 
