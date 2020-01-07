#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM172BP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM172BP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM173A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM173B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM174A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM174B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM177A1 
