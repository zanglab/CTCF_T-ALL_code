#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EPS8L3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EPS8L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPSTI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPYC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EQTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERAL1 
