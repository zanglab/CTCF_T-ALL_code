#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AKAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP14 
