#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TYW3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TYW3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TYW5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i U2AF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i U2AF1L4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i U2AF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i U2SURP 
