#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CA5B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CA5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA5BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA9 
