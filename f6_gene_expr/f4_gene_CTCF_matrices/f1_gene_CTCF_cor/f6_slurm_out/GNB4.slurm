#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GNB4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GNB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNG10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNG11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNG12 
