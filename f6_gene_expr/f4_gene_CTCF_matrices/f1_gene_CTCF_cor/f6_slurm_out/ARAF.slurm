#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARAF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARAF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARCN1 
