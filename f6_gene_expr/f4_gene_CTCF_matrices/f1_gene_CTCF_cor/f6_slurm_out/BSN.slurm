#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BSN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BSN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BSND 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BSPRY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BST2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BSX 
