#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TLE1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TLE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLE4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLE6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLK1 
