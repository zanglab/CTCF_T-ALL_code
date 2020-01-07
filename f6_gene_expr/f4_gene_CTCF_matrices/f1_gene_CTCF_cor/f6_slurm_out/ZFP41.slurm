#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZFP41.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP57 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP64 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP69 
