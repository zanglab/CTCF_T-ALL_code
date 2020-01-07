#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HOXC10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXC10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXC11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXC12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXC13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXC5 
