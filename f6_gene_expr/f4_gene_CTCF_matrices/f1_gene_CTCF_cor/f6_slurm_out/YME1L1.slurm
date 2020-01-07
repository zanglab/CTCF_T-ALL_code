#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o YME1L1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i YME1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YOD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YPEL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YPEL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YPEL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YPEL4 
