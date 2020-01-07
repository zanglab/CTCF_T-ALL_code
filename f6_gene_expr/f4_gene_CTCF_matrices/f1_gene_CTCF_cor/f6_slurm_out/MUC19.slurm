#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MUC19.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC3A 
