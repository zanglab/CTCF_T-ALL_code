#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TTC5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC9 
