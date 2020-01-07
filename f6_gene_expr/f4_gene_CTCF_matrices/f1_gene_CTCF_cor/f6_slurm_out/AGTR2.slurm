#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AGTR2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AGTR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGTRAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGXT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGXT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AHCTF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AHCY 
