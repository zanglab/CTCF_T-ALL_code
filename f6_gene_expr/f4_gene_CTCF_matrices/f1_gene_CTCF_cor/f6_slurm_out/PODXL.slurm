#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PODXL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PODXL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PODXL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POF1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POFUT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POFUT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POGK 
