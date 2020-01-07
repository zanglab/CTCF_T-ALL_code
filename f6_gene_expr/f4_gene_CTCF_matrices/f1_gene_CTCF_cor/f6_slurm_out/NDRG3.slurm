#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NDRG3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NDRG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDRG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDST2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDST3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDST4 
