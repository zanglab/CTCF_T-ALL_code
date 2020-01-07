#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LSM3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LSM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LSM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LSM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LSM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LSM7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LSM8 
