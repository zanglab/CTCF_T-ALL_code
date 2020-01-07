#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LSM10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LSM10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LSM11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LSM12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LSM14A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LSM14B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LSM2 
