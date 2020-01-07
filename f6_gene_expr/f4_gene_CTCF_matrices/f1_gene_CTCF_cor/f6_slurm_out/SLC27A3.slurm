#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC27A3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC27A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC27A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC27A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC27A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC28A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC28A2 
