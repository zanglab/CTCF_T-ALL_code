#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NMT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NMT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NMU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NMUR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NMUR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NNAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NNMT 
