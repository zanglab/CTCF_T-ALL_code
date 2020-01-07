#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARMC8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMCX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMCX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMCX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMCX4 
