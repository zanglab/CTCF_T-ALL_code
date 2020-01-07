#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NTN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NTN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTN5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTNG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTNG2 
