#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM109A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM109A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM109B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM110A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM110B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM110C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM110D 
