#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o M6PR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i M6PR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAATS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAB21L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAB21L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAB21L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MACC1 
