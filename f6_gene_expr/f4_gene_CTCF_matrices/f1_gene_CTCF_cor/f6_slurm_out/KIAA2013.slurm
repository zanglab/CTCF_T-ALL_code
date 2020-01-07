#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIAA2013.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA2013 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA2026 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIDINS220 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF13A 
