#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SRSF8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SRSF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRSF9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRXN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SS18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SS18L1 
