#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NUSAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NUSAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUTF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUTM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUTM2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUTM2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUTM2D 
