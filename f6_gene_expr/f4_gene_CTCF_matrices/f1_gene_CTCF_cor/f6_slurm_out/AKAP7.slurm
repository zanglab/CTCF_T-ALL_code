#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AKAP7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP8L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKIRIN1 
