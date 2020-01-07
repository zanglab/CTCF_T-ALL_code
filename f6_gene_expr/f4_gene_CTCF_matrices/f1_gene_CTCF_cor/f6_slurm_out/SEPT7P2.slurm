#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SEPT7P2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SEPT7P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEPT7P9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEPT8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEPT9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERAC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERBP1 
