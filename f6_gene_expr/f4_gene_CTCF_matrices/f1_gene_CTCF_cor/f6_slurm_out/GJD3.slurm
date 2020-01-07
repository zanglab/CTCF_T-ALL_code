#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GJD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GJD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GK3P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GK5 
