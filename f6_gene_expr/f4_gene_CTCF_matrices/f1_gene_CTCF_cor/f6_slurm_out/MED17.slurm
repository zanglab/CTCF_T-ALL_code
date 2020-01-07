#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MED17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MED17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED22 
