#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BTLA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BTLA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTN1A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTN2A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTN2A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTN2A3P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTN3A1 
