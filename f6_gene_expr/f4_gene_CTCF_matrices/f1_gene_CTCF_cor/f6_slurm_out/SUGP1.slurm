#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SUGP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SUGP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUGP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUGT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUGT1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUGT1P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SULF1 
