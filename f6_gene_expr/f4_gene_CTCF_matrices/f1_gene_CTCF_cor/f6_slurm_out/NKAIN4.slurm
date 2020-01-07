#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NKAIN4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NKAIN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKAPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKAPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKD2 
