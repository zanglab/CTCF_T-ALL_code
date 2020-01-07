#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SFRP5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SFRP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFSWAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFT2D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFT2D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFT2D3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFTA2 
