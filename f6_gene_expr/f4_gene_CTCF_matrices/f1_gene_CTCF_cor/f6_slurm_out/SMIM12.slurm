#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMIM12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM18 
