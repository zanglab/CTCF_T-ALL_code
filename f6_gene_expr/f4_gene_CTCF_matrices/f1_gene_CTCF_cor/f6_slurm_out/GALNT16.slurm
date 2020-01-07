#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GALNT16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT5 
