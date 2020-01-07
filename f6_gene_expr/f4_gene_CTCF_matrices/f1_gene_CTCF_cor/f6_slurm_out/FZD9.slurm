#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FZD9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FZD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FZR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i G0S2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i G2E3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i G3BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i G3BP2 
