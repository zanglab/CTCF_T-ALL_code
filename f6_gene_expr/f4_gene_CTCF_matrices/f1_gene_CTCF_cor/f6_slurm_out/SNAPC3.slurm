#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNAPC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAPC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAPC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAPC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAPIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNCA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNCAIP 
