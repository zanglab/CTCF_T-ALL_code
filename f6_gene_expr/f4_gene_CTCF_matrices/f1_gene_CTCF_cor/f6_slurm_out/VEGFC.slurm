#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VEGFC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VEGFC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VENTX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VENTXP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VEPH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VEZF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VEZT 
