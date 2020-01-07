#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CACNA1S.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNA1S 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNA2D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNA2D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNA2D3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNA2D4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNB1 
