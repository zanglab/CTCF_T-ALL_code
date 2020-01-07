#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ALOX12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ALOX12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALOX12B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALOX12P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALOX15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALOX15B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALOX5 
