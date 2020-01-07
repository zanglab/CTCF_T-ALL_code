#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZFP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP36L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP36L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP37 
