#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NEK3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NEK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEK5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEK6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEK7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEK8 
