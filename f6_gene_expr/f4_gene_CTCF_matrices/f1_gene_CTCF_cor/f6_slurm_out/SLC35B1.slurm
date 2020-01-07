#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC35B1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35B3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35B4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35C2 
