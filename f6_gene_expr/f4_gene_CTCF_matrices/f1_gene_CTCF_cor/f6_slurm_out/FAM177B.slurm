#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM177B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM177B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM178B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM180A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM180B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM181A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM181B 
