#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYH2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYH4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYH6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYH7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYH7B 
