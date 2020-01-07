#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IFI27L1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IFI27L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFI27L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFI30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFI35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFI44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFI44L 
