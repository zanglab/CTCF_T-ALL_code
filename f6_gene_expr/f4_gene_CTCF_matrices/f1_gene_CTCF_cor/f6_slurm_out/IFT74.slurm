#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IFT74.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT74 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT80 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT81 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT88 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGBP1P1 
