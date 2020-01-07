#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLHL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL34 
