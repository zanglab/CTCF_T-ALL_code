#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CASP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP8AP2 
