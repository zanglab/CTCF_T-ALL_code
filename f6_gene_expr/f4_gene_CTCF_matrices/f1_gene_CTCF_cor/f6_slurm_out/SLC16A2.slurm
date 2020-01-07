#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC16A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC16A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC16A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC16A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC16A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC16A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC16A7 
