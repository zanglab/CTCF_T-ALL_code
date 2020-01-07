#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EPRS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EPRS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPS15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPS15L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPS8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPS8L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPS8L2 
