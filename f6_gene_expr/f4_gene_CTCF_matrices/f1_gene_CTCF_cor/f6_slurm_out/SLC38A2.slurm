#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC38A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC38A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC38A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC38A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC38A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC38A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC38A7 
