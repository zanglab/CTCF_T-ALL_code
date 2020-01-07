#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GFAP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GFAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GFER 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GFI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GFI1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GFM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GFM2 
