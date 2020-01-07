#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC5A12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC5A12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC5A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC5A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC5A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC5A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC5A6 
