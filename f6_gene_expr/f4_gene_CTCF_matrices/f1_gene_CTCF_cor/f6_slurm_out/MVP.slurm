#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MVP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MVP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MXD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MXD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MXD4 
