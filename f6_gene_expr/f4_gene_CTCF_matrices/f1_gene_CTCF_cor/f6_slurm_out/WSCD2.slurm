#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WSCD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WSCD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WTAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WTAPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WTIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WWC1 
