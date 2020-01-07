#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABI3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABI3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABI3BP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABLIM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABLIM2 
