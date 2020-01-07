#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HDAC9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HDAC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDDC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDGF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDGFL1 
