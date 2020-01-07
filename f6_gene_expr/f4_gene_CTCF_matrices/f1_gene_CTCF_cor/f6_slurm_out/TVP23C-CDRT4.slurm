#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TVP23C-CDRT4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TVP23C-CDRT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TWF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TWF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TWIST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TWIST2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TWISTNB 
