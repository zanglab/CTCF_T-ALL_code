#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MSC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MSC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSGN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSH4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSH5 
