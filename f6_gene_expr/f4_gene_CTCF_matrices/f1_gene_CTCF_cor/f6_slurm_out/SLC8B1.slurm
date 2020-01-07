#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC8B1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC8B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9A3R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9A3R2 
