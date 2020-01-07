#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FILIP1L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FILIP1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FIP1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FIS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FITM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FITM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FIZ1 
