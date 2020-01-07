#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C1QTNF6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QTNF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QTNF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QTNF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QTNF9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QTNF9B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1R 
