#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLEC4G.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC4G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC4GP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC4M 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC7A 
