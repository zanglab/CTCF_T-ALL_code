#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC38A8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC38A8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC38A9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A12 
