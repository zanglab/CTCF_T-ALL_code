#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PREP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PREP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PREPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PREX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PREX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRG2 
