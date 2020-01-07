#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLK4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLK5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLLP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLOD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLOD2 
