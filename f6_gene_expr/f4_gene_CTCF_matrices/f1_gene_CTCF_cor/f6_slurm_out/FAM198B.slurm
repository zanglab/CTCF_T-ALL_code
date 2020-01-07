#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM198B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM198B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM199X 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM19A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM19A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM19A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM19A4 
