#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZIM2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZIM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZIM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZKSCAN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZKSCAN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZKSCAN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZKSCAN4 
