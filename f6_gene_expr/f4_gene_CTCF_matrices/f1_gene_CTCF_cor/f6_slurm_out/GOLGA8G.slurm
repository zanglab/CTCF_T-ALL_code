#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GOLGA8G.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8J 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8M 
