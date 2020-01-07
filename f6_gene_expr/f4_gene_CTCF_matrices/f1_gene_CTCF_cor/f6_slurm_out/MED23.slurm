#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MED23.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MED23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED28 
