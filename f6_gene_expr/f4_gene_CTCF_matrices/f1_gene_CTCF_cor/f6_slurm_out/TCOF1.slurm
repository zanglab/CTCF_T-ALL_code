#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TCOF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TCOF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCP10L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCP10L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCP11 
