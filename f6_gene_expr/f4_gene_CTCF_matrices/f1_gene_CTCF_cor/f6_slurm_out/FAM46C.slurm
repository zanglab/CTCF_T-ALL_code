#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM46C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM46C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM46D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM47A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM47C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM47E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM47E-STBD1 
