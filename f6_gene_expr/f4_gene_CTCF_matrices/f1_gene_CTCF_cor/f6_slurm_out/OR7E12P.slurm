#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR7E12P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7E12P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7E14P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7E156P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7E24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7E2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7E37P 
