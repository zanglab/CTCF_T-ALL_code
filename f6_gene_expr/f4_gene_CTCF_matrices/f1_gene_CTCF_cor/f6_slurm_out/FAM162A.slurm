#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM162A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM162A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM162B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM163A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM163B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM166A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM166B 
