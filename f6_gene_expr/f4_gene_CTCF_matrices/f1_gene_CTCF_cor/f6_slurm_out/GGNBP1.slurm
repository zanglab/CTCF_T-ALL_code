#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GGNBP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GGNBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGNBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGPS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGT3P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGT5 
