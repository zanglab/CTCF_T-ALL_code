#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC26A10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC26A10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC26A11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC26A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC26A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC26A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC26A5 
