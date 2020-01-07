#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ICAM2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ICAM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ICAM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ICAM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ICAM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ICE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ICE2 
