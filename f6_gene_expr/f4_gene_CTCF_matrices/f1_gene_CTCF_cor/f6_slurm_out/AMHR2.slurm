#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AMHR2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AMHR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMIGO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMIGO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMIGO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMMECR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMMECR1L 
