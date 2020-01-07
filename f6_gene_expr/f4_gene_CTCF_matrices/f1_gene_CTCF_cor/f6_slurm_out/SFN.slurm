#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SFN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SFN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFPQ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFRP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFRP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFRP4 
