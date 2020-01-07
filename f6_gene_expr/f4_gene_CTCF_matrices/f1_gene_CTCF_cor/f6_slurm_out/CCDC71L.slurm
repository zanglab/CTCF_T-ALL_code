#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC71L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC71L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC73 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC74A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC74B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC77 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC78 
