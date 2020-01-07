#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GCSH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GCSH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCSHP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDAP1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDAP2 
