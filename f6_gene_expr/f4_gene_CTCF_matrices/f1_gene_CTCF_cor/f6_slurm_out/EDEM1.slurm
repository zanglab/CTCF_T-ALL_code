#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EDEM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EDEM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDEM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDEM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDIL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDN1 
