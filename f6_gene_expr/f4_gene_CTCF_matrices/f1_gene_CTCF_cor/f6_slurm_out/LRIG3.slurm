#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRIG3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRIG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRIT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRIT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRIT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRP1 
