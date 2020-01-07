#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IDH3B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IDH3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IDH3G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IDI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IDI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IDNK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IDO1 
