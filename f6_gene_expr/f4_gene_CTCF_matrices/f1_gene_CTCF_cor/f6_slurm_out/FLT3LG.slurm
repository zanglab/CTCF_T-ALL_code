#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FLT3LG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FLT3LG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLVCR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLVCR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLYWCH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLYWCH2 
