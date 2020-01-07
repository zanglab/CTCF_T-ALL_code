#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM35A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM35A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM35BP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM35DP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM3C 
