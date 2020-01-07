#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBOX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBQLN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBQLN2 
