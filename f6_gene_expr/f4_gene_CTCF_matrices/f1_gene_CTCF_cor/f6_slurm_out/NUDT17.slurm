#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NUDT17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT22 
