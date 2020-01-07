#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C1orf141.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf141 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf146 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf158 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf159 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf162 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf167 
