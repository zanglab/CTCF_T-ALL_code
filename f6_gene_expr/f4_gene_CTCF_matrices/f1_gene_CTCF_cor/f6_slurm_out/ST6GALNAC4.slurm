#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ST6GALNAC4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ST6GALNAC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST6GALNAC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST6GALNAC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST7L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST8SIA1 
