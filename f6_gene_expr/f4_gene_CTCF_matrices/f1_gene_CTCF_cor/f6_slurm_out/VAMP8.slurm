#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VAMP8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VAMP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VANGL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VANGL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VAPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VAPB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VARS 
