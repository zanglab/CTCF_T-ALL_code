#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C5orf67.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf67 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf106 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf118 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf120 
