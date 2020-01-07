#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAPK3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPK6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPK7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPK8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPK8IP1 
