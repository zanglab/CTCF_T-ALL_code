#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VPS45.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS51 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS52 
