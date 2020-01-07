#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR8D4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR8D4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR8G1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR8G5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR8H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR8H2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR8J3 
