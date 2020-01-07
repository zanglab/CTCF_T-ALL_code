#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TIMP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TINAG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TINAGL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TINCR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TINF2 
