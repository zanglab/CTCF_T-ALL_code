#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM131C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM131C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM133A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM133B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM133CP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM133DP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM135A 
