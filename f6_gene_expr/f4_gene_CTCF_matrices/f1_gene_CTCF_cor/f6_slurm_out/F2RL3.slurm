#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o F2RL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i F2RL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i F3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i F5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i F7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i F8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i F8A1 
