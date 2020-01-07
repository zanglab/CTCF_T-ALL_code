#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IGSF3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IGSF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGSF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGSF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGSF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGSF9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGSF9B 
