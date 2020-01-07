#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C1orf50.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf56 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf61 
