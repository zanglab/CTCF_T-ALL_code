#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ONECUT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ONECUT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ONECUT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OOEP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OOSP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OPA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OPA3 
