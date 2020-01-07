#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HOPX.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HOPX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HORMAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HORMAD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXA10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXA11 
