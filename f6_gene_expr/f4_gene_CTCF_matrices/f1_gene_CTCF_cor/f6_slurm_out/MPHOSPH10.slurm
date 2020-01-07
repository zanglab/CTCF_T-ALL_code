#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MPHOSPH10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MPHOSPH10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPHOSPH6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPHOSPH8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPHOSPH9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPL 
