#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C9orf129.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf129 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf131 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf135 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf152 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf153 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf16 
