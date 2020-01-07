#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLIP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLIP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLIP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLK2P1 
