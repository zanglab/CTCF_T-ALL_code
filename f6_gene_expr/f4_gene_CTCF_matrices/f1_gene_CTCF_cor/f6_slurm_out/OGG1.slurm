#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OGG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OGG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OGN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OGT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OIP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OIT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OLA1 
