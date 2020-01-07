#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYO1C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO1D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO1E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO1F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO1G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO1H 
