#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIFAP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIFAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIFC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIFC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIFC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIR2DL1 
