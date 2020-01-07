#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC28A3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC28A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC29A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC29A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC29A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC29A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC2A1 
