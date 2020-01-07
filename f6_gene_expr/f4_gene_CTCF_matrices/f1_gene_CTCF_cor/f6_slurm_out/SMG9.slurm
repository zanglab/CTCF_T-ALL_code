#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMG9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMG9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM10L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM10L2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM10L2B 
