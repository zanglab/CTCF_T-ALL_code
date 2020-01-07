#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNHD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNLZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNM1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNM1P46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNM2 
