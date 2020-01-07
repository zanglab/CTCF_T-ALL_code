#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MSTO1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MSTO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSTO2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSX2P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MT1A 
