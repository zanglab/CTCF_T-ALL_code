#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIST1H3G.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H3G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H3H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H3I 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H3J 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H4B 
