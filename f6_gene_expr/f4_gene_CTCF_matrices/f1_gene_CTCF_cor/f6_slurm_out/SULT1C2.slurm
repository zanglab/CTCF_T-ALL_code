#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SULT1C2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT1C2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT1C2P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT1C3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT1C4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT1E1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT2A1 
