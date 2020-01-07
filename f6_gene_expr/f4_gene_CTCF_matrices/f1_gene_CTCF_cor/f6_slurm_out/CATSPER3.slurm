#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CATSPER3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CATSPER3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CATSPER4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CATSPERB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CATSPERD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CATSPERG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAV1 
