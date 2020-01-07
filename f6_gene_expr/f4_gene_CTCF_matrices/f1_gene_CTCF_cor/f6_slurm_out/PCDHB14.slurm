#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCDHB14.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB17P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB18P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB19P 
