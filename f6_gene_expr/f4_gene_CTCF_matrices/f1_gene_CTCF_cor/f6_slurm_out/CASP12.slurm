#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CASP12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP16P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP1P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP3 
