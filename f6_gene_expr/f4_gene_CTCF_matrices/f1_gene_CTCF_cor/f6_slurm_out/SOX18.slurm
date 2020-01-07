#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SOX18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX4 
