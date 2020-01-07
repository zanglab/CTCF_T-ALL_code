#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRDM10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM15 
