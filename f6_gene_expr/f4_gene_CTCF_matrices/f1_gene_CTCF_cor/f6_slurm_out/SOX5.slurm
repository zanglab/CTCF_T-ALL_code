#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SOX5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SP1 
