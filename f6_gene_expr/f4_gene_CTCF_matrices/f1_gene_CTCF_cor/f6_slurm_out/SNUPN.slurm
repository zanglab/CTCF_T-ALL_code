#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNUPN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNUPN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNURF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNW1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX11 
