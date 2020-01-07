#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DRC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DRC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRD4 
