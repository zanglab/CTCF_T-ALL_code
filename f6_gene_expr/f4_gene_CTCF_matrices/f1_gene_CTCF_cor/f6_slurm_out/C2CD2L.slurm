#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C2CD2L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C2CD2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2CD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2CD4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2CD4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2CD4C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2CD4D 
