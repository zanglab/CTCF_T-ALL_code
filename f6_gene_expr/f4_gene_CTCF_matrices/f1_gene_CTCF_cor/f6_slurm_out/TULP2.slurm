#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TULP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TULP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TULP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TULP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUNAR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUSC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUSC2 
