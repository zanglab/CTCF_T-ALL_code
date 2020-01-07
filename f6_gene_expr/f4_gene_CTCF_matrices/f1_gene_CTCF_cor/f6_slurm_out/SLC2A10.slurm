#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC2A10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC2A10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC2A11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC2A12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC2A13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC2A14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC2A2 
