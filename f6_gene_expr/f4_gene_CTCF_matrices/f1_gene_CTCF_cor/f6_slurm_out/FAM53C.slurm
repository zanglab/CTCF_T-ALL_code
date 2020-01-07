#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM53C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM53C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM57A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM57B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM69A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM69B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM69C 
