#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FLOT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FLOT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLRT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLRT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLRT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLT3 
