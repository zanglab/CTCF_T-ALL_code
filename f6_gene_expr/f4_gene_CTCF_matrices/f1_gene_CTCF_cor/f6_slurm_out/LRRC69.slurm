#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRRC69.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC69 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC70 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC71 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC72 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC73 
