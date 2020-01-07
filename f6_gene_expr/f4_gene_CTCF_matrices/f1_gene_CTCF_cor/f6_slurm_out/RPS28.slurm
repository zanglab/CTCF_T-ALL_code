#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPS28.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS2P32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS4X 
