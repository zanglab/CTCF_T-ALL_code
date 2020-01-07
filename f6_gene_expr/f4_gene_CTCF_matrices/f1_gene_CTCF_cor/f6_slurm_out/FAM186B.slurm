#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM186B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM186B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM187B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM189A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM189A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM189B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM192A 
