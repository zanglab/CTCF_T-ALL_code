#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GUSBP10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GUSBP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUSBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUSBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUSBP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUSBP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUSBP9 
