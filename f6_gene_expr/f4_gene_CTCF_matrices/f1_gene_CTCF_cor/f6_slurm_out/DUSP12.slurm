#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DUSP12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP18 
