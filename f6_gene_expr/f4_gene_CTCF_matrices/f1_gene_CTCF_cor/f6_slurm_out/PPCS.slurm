#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPCS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPCS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPDPF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPEF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPEF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPFIA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPFIA2 
