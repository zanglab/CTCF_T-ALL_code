#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM160A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM160A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM160A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM160B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM160B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM161A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM161B 
