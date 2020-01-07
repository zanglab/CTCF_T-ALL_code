#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C1QTNF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QTNF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QTNF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QTNF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QTNF3-AMACR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QTNF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QTNF5 
