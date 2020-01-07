#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TBC1D8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D8B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D9B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBCA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBCB 
