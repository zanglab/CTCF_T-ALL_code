#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MED29.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MED29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED7 
