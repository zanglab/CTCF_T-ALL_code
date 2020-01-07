#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRPS1L1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPS1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPSAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPSAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR12 
