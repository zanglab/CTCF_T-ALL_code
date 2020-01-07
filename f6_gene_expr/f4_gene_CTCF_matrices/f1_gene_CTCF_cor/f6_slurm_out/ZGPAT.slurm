#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZGPAT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZGPAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZGRF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZHX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZHX1-C8orf76 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZHX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZHX3 
