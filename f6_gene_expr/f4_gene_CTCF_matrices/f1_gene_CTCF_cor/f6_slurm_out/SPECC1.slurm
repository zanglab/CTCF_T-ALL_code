#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPECC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPECC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPECC1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPECC1L-ADORA2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPEF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPEF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPEG 
