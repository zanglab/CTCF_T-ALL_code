#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM193A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM193A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM193B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM196A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM196B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM197Y2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM198A 
