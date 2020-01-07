#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZSWIM7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSWIM7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSWIM8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZUFSP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZW10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZWILCH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZWINT 
