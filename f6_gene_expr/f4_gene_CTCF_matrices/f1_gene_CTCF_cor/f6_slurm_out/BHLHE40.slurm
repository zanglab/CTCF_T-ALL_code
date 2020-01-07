#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BHLHE40.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BHLHE40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BHLHE41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BHMG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BHMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BHMT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BICC1 
