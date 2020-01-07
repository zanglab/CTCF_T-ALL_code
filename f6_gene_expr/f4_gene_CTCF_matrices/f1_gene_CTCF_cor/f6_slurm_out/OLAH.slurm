#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OLAH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OLAH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OLFM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OLFM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OLFM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OLFM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OLFML1 
