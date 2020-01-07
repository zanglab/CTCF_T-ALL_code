#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DHFR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DHFR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHODH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHPS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHRS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHRS11 
