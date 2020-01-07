#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC2A3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC2A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC2A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC2A4RG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC2A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC2A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC2A7 
