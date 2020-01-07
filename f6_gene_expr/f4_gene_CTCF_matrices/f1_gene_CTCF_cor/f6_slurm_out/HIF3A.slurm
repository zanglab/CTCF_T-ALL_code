#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIF3A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIF3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIGD1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIGD1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIGD1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIGD2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIGD2B 
