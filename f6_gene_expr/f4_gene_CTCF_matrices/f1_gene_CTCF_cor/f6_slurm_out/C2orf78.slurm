#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C2orf78.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf78 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf80 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf81 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf83 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf88 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf91 
