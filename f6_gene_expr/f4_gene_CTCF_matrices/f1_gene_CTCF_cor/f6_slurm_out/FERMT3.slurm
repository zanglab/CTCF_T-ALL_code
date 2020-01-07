#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FERMT3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FERMT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FES 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FETUB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FEV 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FEZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FEZ2 
