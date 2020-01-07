#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OS9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OS9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSBPL10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSBPL11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSBPL1A 
