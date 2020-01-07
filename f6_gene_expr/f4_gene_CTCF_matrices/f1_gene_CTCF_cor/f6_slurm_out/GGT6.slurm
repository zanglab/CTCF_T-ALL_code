#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GGT6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GGT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGT7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGT8P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGTA1P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGTLC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGTLC2 
