#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC4A11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC4A11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC4A1AP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC4A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC4A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC4A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC4A5 
