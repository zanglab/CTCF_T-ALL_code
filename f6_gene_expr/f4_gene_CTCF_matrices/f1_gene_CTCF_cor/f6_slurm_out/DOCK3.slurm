#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DOCK3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DOCK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOCK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOCK5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOCK6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOCK7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOCK8 
