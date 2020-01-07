#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C4B_2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C4B_2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4orf17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4orf19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4orf22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4orf3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4orf33 
