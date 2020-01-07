#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CTDSPL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CTDSPL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTGF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTHRC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTIF 
