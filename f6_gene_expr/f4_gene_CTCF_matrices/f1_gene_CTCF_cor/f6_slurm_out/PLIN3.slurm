#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLIN3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLIN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLIN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLIN5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLK3 
