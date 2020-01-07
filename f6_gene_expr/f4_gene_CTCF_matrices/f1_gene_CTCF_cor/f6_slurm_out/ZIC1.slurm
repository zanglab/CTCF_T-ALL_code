#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZIC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZIC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZIC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZIC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZIC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZIC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZIK1 
