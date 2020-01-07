#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VSIG2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VSIG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VSIG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VSIG8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VSNL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VSTM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VSTM2A 
