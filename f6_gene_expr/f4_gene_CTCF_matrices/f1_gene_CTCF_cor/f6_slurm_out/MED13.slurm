#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MED13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MED13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED13L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED15P9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED16 
