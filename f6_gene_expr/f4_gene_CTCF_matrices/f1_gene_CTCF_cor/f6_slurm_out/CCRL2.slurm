#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCRL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCRL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCSAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCSER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCSER2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCT2 
