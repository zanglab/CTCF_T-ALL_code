#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTK7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTK7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTMA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTMS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTOV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTP4A1 
