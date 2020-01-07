#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM71F1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM71F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM71F2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM72A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM72B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM72C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM72D 
