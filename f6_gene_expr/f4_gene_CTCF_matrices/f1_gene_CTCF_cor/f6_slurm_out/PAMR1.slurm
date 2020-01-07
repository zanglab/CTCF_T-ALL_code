#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PAMR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PAMR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PANK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PANK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PANK3 
