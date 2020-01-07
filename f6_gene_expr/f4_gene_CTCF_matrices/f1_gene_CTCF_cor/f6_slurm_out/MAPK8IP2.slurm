#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAPK8IP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPK8IP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPK8IP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPK9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPKAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPKAPK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPKAPK3 
