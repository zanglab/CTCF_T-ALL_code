#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MIPOL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MIPOL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIS12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIS18A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIS18BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MISP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MITD1 
