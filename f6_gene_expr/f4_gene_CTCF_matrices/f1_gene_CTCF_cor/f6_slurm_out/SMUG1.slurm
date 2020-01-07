#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMUG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMUG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMURF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMURF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMYD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMYD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMYD3 
