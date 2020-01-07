#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC9A9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9A9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC9C2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLCO1A2 
