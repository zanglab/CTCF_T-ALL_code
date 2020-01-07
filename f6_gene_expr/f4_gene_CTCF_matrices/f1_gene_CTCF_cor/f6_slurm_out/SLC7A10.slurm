#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC7A10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A3 
