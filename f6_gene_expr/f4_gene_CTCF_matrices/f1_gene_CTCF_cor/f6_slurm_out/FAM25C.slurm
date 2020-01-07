#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM25C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM25C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM25G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM27B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM27E2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM27E3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM32A 
