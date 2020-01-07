#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RBMS3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RBMS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBMX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBMX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBMXL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBMXL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBMXL3 
