#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLEC1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC3A 
