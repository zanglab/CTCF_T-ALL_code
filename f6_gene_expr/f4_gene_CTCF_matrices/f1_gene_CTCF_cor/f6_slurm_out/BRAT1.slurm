#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BRAT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BRAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRCA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRCA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRCC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRD2 
