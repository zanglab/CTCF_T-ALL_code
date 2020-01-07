#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC37A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC37A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC37A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC37A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC38A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC38A10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC38A11 
