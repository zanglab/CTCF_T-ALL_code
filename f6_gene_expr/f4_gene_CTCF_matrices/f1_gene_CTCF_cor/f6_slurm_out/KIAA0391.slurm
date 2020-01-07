#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIAA0391.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0391 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0408 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0513 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0556 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0586 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0753 
