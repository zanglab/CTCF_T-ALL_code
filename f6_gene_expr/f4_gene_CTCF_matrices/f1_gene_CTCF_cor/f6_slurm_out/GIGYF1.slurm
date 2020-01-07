#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GIGYF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GIGYF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GIGYF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GIMAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GIMAP1-GIMAP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GIMAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GIMAP4 
