#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM103A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM103A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM104A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM104B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM105A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM107A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM107B 
