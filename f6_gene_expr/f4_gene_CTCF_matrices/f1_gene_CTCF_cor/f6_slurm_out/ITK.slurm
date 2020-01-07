#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ITK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ITK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITLN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITLN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITM2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITM2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITM2C 
