#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIAA0100.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0100 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0141 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0232 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0319 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0319L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0355 
