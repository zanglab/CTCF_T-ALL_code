#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SYDE1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SYDE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYDE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYMPK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYN1 
