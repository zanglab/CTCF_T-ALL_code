#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SIKE1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SIKE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIMC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIN3A 
