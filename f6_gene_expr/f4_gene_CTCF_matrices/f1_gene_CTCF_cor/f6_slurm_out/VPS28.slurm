#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VPS28.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS33A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS33B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS36 
