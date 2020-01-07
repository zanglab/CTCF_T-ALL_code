#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RGS3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS7BP 
