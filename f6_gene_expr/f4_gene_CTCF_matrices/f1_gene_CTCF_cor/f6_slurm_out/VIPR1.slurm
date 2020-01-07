#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VIPR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VIPR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VIPR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VIT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VKORC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VKORC1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VLDLR 
