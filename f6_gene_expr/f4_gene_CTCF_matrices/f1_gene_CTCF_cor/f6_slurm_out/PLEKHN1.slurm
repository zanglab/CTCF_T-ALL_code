#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLEKHN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLET1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLG 
