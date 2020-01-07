#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAST1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAST2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAST3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAST4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MASTL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAT1A 
