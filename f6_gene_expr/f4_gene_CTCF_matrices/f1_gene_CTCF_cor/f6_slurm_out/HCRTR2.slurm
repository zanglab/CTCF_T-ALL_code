#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HCRTR2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HCRTR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HCST 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDAC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDAC10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDAC11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDAC2 
