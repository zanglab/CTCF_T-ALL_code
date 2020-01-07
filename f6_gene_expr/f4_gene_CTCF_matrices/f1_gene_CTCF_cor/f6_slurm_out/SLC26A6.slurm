#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC26A6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC26A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC26A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC26A8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC26A9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC27A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC27A2 
