#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARMC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMC7 
