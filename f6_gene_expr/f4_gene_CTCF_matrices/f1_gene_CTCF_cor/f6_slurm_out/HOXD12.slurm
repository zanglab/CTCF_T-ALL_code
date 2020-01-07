#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HOXD12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXD12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXD13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXD9 
