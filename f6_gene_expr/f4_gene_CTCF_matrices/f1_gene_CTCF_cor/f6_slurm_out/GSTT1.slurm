#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GSTT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTT2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSX2 
