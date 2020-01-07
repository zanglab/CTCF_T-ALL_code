#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IFT22.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT57 
