#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZBTB7C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB7C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB8A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB8B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB8OS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC2HC1A 
