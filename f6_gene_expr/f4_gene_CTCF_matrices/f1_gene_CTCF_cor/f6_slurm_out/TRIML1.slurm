#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRIML1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIML1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIML2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIOBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIP11 
