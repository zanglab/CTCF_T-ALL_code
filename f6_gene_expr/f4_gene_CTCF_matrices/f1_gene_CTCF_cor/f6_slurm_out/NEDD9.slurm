#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NEDD9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NEDD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEFH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEFL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEFM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEGR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEIL1 
