#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o B4GALT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i B4GALT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B4GALT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B4GALT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B4GALT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B4GALT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B4GALT7 
