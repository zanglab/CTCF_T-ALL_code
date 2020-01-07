#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EMX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EMX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENAH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENAM 
