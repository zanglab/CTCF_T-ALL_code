#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IFNA4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNAR2 
