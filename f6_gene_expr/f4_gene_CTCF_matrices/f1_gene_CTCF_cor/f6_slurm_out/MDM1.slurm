#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MDM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MDM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MDM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MDM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MDN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MDP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ME1 
