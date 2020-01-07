#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GRPR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GRPR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRSF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRTP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRWD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRXCR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRXCR2 
