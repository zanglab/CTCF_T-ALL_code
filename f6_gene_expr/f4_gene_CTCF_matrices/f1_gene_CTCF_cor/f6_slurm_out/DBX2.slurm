#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DBX2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DBX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCAF10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCAF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCAF12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCAF12L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCAF12L2 
