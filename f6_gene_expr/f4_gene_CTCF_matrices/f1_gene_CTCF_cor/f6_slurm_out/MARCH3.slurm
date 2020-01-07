#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MARCH3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MARCH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARCH4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARCH5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARCH6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARCH7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARCH8 
