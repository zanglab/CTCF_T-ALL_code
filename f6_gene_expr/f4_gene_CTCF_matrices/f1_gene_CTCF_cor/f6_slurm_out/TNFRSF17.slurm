#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TNFRSF17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF21 
