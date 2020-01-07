#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NFATC4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NFATC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFE2L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFE2L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFE2L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFE4 
