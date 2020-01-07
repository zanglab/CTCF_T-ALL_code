#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABTB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABTB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABTB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACAA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACAA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACACA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACACB 
