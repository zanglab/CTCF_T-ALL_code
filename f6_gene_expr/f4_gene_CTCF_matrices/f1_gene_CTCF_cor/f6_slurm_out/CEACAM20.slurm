#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CEACAM20.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CEACAM20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEACAM21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEACAM22P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEACAM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEACAM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEACAM5 
