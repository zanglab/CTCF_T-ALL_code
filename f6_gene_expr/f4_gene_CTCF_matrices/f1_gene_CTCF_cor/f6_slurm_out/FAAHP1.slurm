#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAAHP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAAHP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAAP100 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAAP20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAAP24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FABP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FABP12 
