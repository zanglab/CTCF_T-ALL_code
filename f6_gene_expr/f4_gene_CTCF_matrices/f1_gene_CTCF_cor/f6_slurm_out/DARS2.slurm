#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DARS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DARS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAW1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAXX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAZ2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAZ3 
