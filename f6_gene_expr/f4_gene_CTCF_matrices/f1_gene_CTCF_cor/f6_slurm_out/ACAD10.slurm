#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACAD10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACAD10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACAD11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACAD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACAD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACADL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACADM 
