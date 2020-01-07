#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IL1RAPL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IL1RAPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL1RAPL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL1RL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL1RL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL1RN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL2 
