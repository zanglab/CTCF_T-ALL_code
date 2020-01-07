#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COMMD10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COMMD10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COMMD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COMMD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COMMD3-BMI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COMMD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COMMD5 
