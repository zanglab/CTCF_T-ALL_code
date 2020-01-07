#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM71A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM71A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM71B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM71C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM71D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM71E1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM71E2 
