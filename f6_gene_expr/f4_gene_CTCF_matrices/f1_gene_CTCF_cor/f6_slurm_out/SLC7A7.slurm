#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC7A7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC8A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC8A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC8A3 
