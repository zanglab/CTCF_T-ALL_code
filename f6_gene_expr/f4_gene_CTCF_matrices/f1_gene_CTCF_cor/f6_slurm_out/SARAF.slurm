#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SARAF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SARAF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SARDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SARM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SARNP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SARS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SARS2 
