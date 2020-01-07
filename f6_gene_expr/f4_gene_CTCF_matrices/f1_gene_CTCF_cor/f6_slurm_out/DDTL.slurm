#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DDTL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DDTL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX11L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX11L10 
