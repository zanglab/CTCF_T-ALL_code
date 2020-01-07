#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C6orf132.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf132 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf136 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf141 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf163 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf201 
