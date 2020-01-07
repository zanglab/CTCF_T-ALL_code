#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM167A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM167A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM167B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM168A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM168B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM169A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM169B 
