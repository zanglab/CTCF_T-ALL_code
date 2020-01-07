#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM49A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM49A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM49B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM50A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM50B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM53A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM53B 
