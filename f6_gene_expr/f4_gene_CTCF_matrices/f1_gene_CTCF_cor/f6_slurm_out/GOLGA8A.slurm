#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GOLGA8A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8CP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8DP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8EP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8F 
