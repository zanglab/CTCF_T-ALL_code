#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ERLIN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ERLIN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERMAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERMARD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERMN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERN1 
