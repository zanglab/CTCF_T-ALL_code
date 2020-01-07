#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VDAC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VDAC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VDAC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VDAC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VDR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VEGFA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VEGFB 
