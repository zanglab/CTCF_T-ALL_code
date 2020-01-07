#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POTEA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POTEB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POTEB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POTEC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POTED 
