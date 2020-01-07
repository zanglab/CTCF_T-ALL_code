#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CTCFL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CTCFL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTDNEP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTDP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTDSP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTDSP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTDSPL 
