#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRRC74A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC74A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC74B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC75A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC75B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC8A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC8B 
