#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GCM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GCM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCNT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCNT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCNT3 
