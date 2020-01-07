#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MT1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MT1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MT1DP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MT1E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MT1F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MT1G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MT1H 
