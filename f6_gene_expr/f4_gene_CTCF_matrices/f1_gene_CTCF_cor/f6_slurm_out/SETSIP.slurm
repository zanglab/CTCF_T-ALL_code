#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SETSIP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SETSIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SETX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEZ6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEZ6L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEZ6L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SF1 
